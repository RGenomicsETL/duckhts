/*
 * duckvep_workspace_internal.h — private workspace inspection/plumbing.
 *
 * Not part of the stable kernel ABI. Tests and in-kernel production routing use
 * this to borrow workspace-owned scratch without exposing the
 * opaque workspace layout through kernel/include/duckvep_kernel.h. The returned
 * scratch is sized for one small-variant CodingContext build under the uint16_t allele
 * payload model.
 * Route stats are opt-in internal evidence for tests and are disabled until reset.
 */
#ifndef DUCKVEP_WORKSPACE_INTERNAL_H
#define DUCKVEP_WORKSPACE_INTERNAL_H

#include "duckvep_delta.h"
#include "duckvep_kernel.h"

typedef struct duckvep_workspace_delta_route_stats {
    uint64_t substitution_context;
    uint64_t del_context;
    uint64_t ins_context;
    uint64_t indel_context;
    uint64_t boundary_context;
} duckvep_workspace_delta_route_stats_t;

DUCKVEP_INTERNAL_API duckvep_delta_scratch_t *duckvep_workspace_delta_scratch(
    duckvep_workspace_t *workspace);

DUCKVEP_INTERNAL_API void duckvep_workspace_delta_route_stats_reset(
    duckvep_workspace_t *workspace);
DUCKVEP_INTERNAL_API const duckvep_workspace_delta_route_stats_t *
duckvep_workspace_delta_route_stats(const duckvep_workspace_t *workspace);

#endif /* DUCKVEP_WORKSPACE_INTERNAL_H */
