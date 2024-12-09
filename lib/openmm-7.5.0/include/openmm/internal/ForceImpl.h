#ifndef OPENMM_FORCEIMPL_H_
#define OPENMM_FORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/internal/windowsExport.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace OpenMM {

class Force;
class ContextImpl;

/**
 * A ForceImpl provides the internal implementation of a Force.  When a Context is
 * created for a System, it creates a ForceImpl for each Force in the System.  The ForceImpl
 * is permitted to cache information from the Force or System which is needed to apply the
 * force efficiently.  If the user calls reinitialize() on the Context, all ForceImpl
 * objects it had previously created are deleted and recreated from scratch.
 * 
 * This is an abstract class.  Each Force subclass is responsible for defining its own
 * ForceImpl subclass.
 */

class OPENMM_EXPORT ForceImpl {
public:
    virtual ~ForceImpl() {
    }
    /**
     * This is called after the ForceImpl is created and before updateContextState(), calcForces(),
     * or calcEnergy() is called on it.  This allows it to do any necessary initialization.
     */
    virtual void initialize(ContextImpl& context) = 0;
    /**
     * Get the Force object from which this ForceImpl was created.
     */
    virtual const Force& getOwner() const = 0;
    /**
     * This method is called at the beginning of each time step.  It give the ForceImpl a chance
     * to modify the state variables (positions, velocities, and parameters) stored in the
     * Context in arbitrary ways before integration is performed.
     * 
     * @param context        the context in which the system is being simulated
     * @param forcesInvalid  if the state was modified in any way that might cause previously
     *                       calculated forces to no longer be valid (such as modifying
     *                       positions or parameters), the method should set this to true.
     */
    virtual void updateContextState(ContextImpl& context, bool& forcesInvalid);
    /**
     * @deprecated This version exists for backward compatibility.  Subclasses should implement the other version instead.
     */
    virtual void updateContextState(ContextImpl& context);
    /**
     * Calculate the force on each particle generated by this ForceImpl and/or this ForceImpl's
     * contribution to the potential energy of the system.
     * 
     * @param context        the context in which the system is being simulated
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param groups         a set of bit flags for which force groups to include.  
     *                       Group i in state j will be included if (groups[j]&(1<<(i-32*j)) != 0.
     * @return this force's contribution to the potential energy of the system, or 0 if this
     * force does not contribute to potential energy (or if includeEnergy is false)
     */
    virtual double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, const std::vector<int>& groups) = 0;
    /**
     * Get a map containing the default values for all adjustable parameters defined by this ForceImpl.  These
     * parameters and their default values will automatically be added to the Context.
     */
    virtual std::map<std::string, double> getDefaultParameters() = 0;
    /**
     * Get the names of all Kernels used by this Force.
     */
    virtual std::vector<std::string> getKernelNames() = 0;
    /**
     * Get pairs of particles connected by bonds by this force.  This is used to determine which particles
     * are part of the same molecule.
     */
    virtual std::vector<std::pair<int, int> > getBondedParticles() const {
        return std::vector<std::pair<int, int> >(0);
    }
protected:
    /**
     * Get the ContextImpl corresponding to a Context.
     */
    ContextImpl& getContextImpl(Context& context) {
        return context.getImpl();
    }
};

} // namespace OpenMM

#endif /*OPENMM_FORCEIMPL_H_*/
