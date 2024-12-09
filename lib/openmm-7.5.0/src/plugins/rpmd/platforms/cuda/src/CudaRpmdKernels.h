#ifndef CUDA_RPMD_KERNELS_H_
#define CUDA_RPMD_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2018 Stanford University and the Authors.      *
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

#include "openmm/RpmdKernels.h"
#include "CudaContext.h"
#include "CudaArray.h"
#include <map>

namespace OpenMM {

/**
 * This kernel is invoked by RPMDIntegrator to take one time step, and to get and
 * set the state of system copies.
 */
class CudaIntegrateRPMDStepKernel : public IntegrateRPMDStepKernel {
public:
    CudaIntegrateRPMDStepKernel(const std::string& name, const Platform& platform, CudaContext& cu) :
            IntegrateRPMDStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the RPMDIntegrator this kernel will be used for
     */
    void initialize(const System& system, const RPMDIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated
     */
    void execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid);
    /**
     * Compute the kinetic energy.
     * 
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator);
    /**
     * Get the positions of all particles in one copy of the system.
     */
    void setPositions(int copy, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles in one copy of the system.
     */
    void setVelocities(int copy, const std::vector<Vec3>& velocities);
    /**
     * Copy positions and velocities for one copy into the context.
     */
    void copyToContext(int copy, ContextImpl& context);
private:
    void computeForces(ContextImpl& context);
    std::string createFFT(int size, const std::string& variable, bool forward);
    CudaContext& cu;
    int numCopies, numParticles, workgroupSize;
    //std::map<int, int> groupsByCopies;
    //int groupsNotContracted;
    std::map<int, std::vector<int>> groupsByCopies;
    std::vector<int> groupsNotContracted;
    CudaArray forces;
    CudaArray positions;
    CudaArray velocities;
    CudaArray contractedForces;
    CudaArray contractedPositions;
    CUfunction pileKernel, stepKernel, velocitiesKernel, copyToContextKernel, copyFromContextKernel, translateKernel;
    std::map<int, CUfunction> positionContractionKernels;
    std::map<int, CUfunction> forceContractionKernels;
};

} // namespace OpenMM

#endif /*CUDA_RPMD_KERNELS_H_*/
