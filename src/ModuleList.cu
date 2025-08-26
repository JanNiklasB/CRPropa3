#include "crpropa/ModuleList.h"

namespace crpropa{

__global__ void cudarun(const thrust::device_vector<Candidate*>& candidates, const ModuleList* MLIST,
	bool recursive, bool secondariesFirst) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i>=candidates.size()) return;

	Candidate* candidate = candidates[i];

	// propagate primary candidate until finished
	while (candidate->isActive()) {
		MLIST->process(candidate);

		// propagate all secondaries before next step of primary
		if (recursive and secondariesFirst) {
			int threadsPerBlock = MLIST->maxThreadsPerBlock;
			int blocksPerGrid   = (candidate->secondaries.size() + threadsPerBlock - 1) / threadsPerBlock;
			cudarun<<<blocksPerGrid, threadsPerBlock>>>(candidate->secondaries, MLIST, recursive, secondariesFirst);
		}
	}

	// propagate secondaries after completing primary
	if (recursive and not secondariesFirst) {
		int threadsPerBlock = MLIST->maxThreadsPerBlock;
		int blocksPerGrid   = (candidate->secondaries.size() + threadsPerBlock - 1) / threadsPerBlock;
		cudarun<<<blocksPerGrid, threadsPerBlock>>>(candidate->secondaries, MLIST, recursive, secondariesFirst);
	}

}

void ModuleList::run(SourceInterface *source, size_t count, bool recursive, bool secondariesFirst) {

	// extract all primaries:
	candidate_vector_t candidates;
	for (size_t i = 0; i < count; i++) {
		ref_ptr<Candidate> candidate;

		candidate = source->getCandidate();

		if (candidate.valid()) {
			candidates.push_back(candidate);
		}
	}

	// start all primaries on gpu:
	int threadsPerBlock = maxThreadsPerBlock;
	int blocksPerGrid   = (candidates.size() + threadsPerBlock - 1) / threadsPerBlock;
	cudarun<<<blocksPerGrid, threadsPerBlock>>>(candidates, this, recursive, secondariesFirst);
}

}