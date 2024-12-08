/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsOpenMM.h"

void DynamicsOpenMM::init() {
    DynamicsBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for OpenMM.");
    if (dyn_type.substr(0,6) != "OpenMM")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for OpenMM.");
    // Create smarter pointer to OpenMM intergrator
    ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    integrator = ha->integrator;
    // Get time step size (in ps) from OpenMM integrator.
    DT = integrator->getStepSize();
}

void DynamicsOpenMM::beforeOneTraj() {
    // Load R,V from trajectory if nucl_load is not empty, otherwise do nothing.
    samplingNucl();
    // Get current time (in ps) and step from OpenMM Context
    // This may not zero, if it is a restart simulation.
    // Here, round() will return the nearest integer value.
    // For multi-traj simulation, they are zero.
    step = round(ha->getTime()/DT);
}

void DynamicsOpenMM::dynamics(int steps) {
    if (dyn_type == "OpenMM") // original OpenMM propagtaion
        integrator->step(steps);
    else // propagation with external forces.
        for (int i = 0; i < steps; ++i)
            myOneStep();
    this->step += steps;
}

void DynamicsOpenMM::myOneStep() {
    static const int propagate_state = param.getInt("propagate_state");
    // Use CompoundIntegrator with two CustomIntegrator to do velocity Verlet
    // integrator with external forces.
    if (param.getStr("integrator") == "velocityVerlet") {
        static std::shared_ptr<OpenMM::CompoundIntegrator> compound =
            std::static_pointer_cast<OpenMM::CompoundIntegrator>(integrator);
        // Part 1 of velocity Verlet integrator:
        // update V with first half step and the R with a full step and the
        // distances constraints.
        // We need compute force at initial step or positions are changed by
        // barostat in ha->updateContextState() [return true if changed]
        // If the positions are the same, then we don't need to compute force here.
        if (ha->updateContextState() || step == 0) {
            ha->getPotentialEnergy(propagate_state, true);
            ha->getForces();
        }
        compound->setCurrentIntegrator(0);
        compound->oneStep(ha->F);
        // Part 2 of velocity Verlet integrator:
        // update V with second half step with new forces: v = v+0.5*dt*f/m
        // and the velocities constraints.
        // Calculate forces based on the updated R
        ha->getPotentialEnergy(propagate_state, true);
        ha->getForces();
        compound->setCurrentIntegrator(1);
        compound->oneStep(ha->F);
        // Note: we need to modify the time/step in platform data since in the
        // above strategy, the simulation was propagated with two steps.
        ha->setTime(ha->getTime() - DT);
        ha->setStep(ha->getStep() - 1);
    }
    else { // The original OpenMM integrator
        ha->updateContextState();
        ha->getPotentialEnergy(propagate_state, true);
        ha->getForces();
        integrator->oneStep(ha->F);
    }
}