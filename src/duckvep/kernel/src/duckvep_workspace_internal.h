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
    uint64_t far_directional_snv;
    uint64_t simple_point_snv;
    uint64_t generalized_pair;
    uint64_t simple_indel;
    uint64_t substitution_context;
    uint64_t del_context;
    uint64_t ins_context;
    uint64_t indel_context;
    uint64_t boundary_context;
    uint64_t uploaded_feature_context;
} duckvep_workspace_delta_route_stats_t;

DUCKVEP_INTERNAL_API duckvep_delta_scratch_t *duckvep_workspace_delta_scratch(
    duckvep_workspace_t *workspace);

DUCKVEP_INTERNAL_API void duckvep_workspace_delta_route_stats_reset(
    duckvep_workspace_t *workspace);
DUCKVEP_INTERNAL_API const duckvep_workspace_delta_route_stats_t *
duckvep_workspace_delta_route_stats(const duckvep_workspace_t *workspace);

/* Test-only authority switch.  Production leaves this disabled.  Enabling it
 * keeps candidate discovery unchanged while routing transcript pairs through
 * the exhaustive region/splice path, so optimized sorted-stream output can be
 * compared field-for-field with the generalized implementation. */
DUCKVEP_INTERNAL_API void duckvep_workspace_force_generalized_annotation(
    duckvep_workspace_t *workspace,
    int                   enabled);

#endif /* DUCKVEP_WORKSPACE_INTERNAL_H */
