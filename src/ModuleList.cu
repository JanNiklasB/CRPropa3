#include "ModuleList.h"

namespace crpropa{

void cudarun(const thrust::device_vector<ref_ptr<Candidate>>& candidates, const ref_ptr<ModuleList> ModuleList,
	bool recursive, bool secondariesFirst) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i>=candidates.size()) return;

	Candidate* candidate = thrust::raw_pointer_cast(candidates[i]);

	// propagate primary candidate until finished
	while (candidate->isActive()) {
		ModuleList->process(candidate);

		// propagate all secondaries before next step of primary
		if (recursive and secondariesFirst) {
			thrust::device_vector<ref_ptr<Candidate>> secondaries = candidate->secondaries;
			int threadsPerBlock = ModuleList->maxThreadsPerBlock;
			int blocksPerGrid   = (secondaries.size() + threadsPerBlock - 1) / threadsPerBlock;
			cudarun<<<blocksPerGrid, threadsPerBlock>>>(secondaries, recursive, secondariesFirst);
		}
	}

	// propagate secondaries after completing primary
	if (recursive and not secondariesFirst) {
		thrust::device_vector<ref_ptr<Candidate>> secondaries = candidate->secondaries;
		int threadsPerBlock = ModuleList->maxThreadsPerBlock;
		int blocksPerGrid   = (secondaries.size() + threadsPerBlock - 1) / threadsPerBlock;
		cudarun<<<blocksPerGrid, threadsPerBlock>>>(secondaries, recursive, secondariesFirst);
	}

}

void ModuleList::cudarun(SourceInterface *source, size_t count, bool recursive, bool secondariesFirst) {

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
	cudarun<<<blocksPerGrid, threadsPerBlock>>>(candidates, recursive, secondariesFirst);
}

}