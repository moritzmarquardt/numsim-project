#pragma once

#include "computation.hpp"
#include "partitioning/partitioning.hpp"
#include <mpi.h>
#include "output_writer/output_writer_paraview_parallel.hpp"
#include "output_writer/output_writer_text_parallel.hpp"

class ParallelComputation : public Computation {
    public:
        void initialize(int argc, char *argv[]) override;
        void runSimulation() override;

    protected:
        void applyBoundaryValues() override;
        void applyInitialBoundaryValues() override;
        void computeTimeStepWidth() override;

        std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaview_;
        std::unique_ptr<OutputWriterTextParallel> outputWriterText_;
};
