/*
 * Copyright (c) 2017 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Parsing DPGP
/** \file
 * \author Anthony J. Greenberg
 *
 * Extracting SNPs from the DPGP .seq files. The variant table will be in the _plink_ BED format. Each chromosome is processed by its own thread in parallel.
 */

#include "sequence.hpp"
#include <vector>
#include <string>
#include <thread>


using std::vector;
using std::string;
using std::thread;

int main(){
	
	vector<string> chromIDs {"Chr2L", "Chr2R", "Chr3L", "Chr3R", "ChrX"}; // chromosome IDs
	vector<unsigned short> chromNums {2, 3, 4, 5, 1};                     // chromosome numbers (needed for the BED metadata)
	vector<thread> threads(4); // main thread will do the X
	
	auto chrIt  = chromIDs.begin();
	auto chrNit = chromNums.begin();
	auto thrIt  = threads.begin();
	
	for (; thrIt != threads.end(); ++thrIt, ++chrIt, ++chrNit) {
		const string inFlList = "seqList_" + (*chrIt) + ".txt";
		const string outFl    = "snp_" + (*chrIt) + ".bed";
		// parsing the autosomes
		SFparse parseA(inFlList, outFl, *chrIt, *chrNit, "SEQ", "BED");
		
		(*thrIt) = thread(parseA);
	}
	const string inFlList("seqList_ChrX.txt");
	const string outFl("snp_ChrX.bed");
	
	// parsing the X
	SFparse parseX(inFlList, outFl, "ChrX", 1, "SEQ", "BED");
	parseX();
	
	for (auto thrdIt = threads.begin(); thrdIt != threads.end(); ++thrdIt) {
		if (thrdIt->joinable()) {
			thrdIt->join();
		}
	}
	
	
}
