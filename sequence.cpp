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

/// Sequence and SNP file parsing and conversion
/** \file
 * \author Anthony J. Greenberg
 * \version 0.9
 *
 * Implementation of facilities that deal with common sequence and variant table file types.
 *
 */


#include "sequence.hpp"
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>

using std::vector;
using std::unordered_map;
using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::numeric_limits;
using std::ceil;

SFparse::SFparse(const vector<string> &inFlNam, const vector<string> &lineNames, const string &refFlNam, const string &outFlNam, const string &chrNam, const unsigned short &chrNum, const string &inFlType, const string &outFlType, const unsigned long &alloc) : _inFileNames(inFlNam), _lineNames(lineNames), _refFlName(refFlNam), _outFileName(outFlNam), _chromName(chrNam), _chromNum(chrNum), _inFileType(inFlType), _outFileType(outFlType), _bufAlloc(alloc) {
	auto outIt = _outFileName.end();
	if ( *(outIt - 4) == '.') { // there is a potentially valid extension
		_outFileName.erase(outIt - 4, outIt); // erase the extension
	}
}

SFparse::SFparse(const string &fileList, const string &outFlNam, const string &chrNam, const unsigned short &chrNum, const string &inFlType, const string &outFlType, const unsigned long &alloc) : _outFileName(outFlNam), _chromName(chrNam), _chromNum(chrNum), _inFileType(inFlType), _outFileType(outFlType), _bufAlloc(alloc) {
	auto outIt = _outFileName.end();
	if ( *(outIt - 4) == '.') { // there is a potentially valid extension
		_outFileName.erase(outIt - 4, outIt); // erase the extension
	}
	
	ifstream lstIn(fileList.c_str());
	
	if (!lstIn) {
		cerr << "ERROR: cannot open file " << fileList << " listing files to process in SFsparse constructor" << endl;
		exit(1);
	}
	string eachFile;
	while (getline(lstIn, eachFile)) {
		if ( (eachFile[0] == 'r') && (eachFile[1] == ':') ) {
			eachFile.erase(0, 2);
			_refFlName = eachFile;
		} else {
			_inFileNames.push_back(eachFile);
		}
		
	}
	lstIn.close();
	if (_inFileNames.size() > numeric_limits<unsigned short>::max()) {
		cerr << "WARNING: number of lines " << _inFileNames.size() << " larger than allowed (" << numeric_limits<unsigned short>::max() << ")" << endl;
	}
	
	for (auto flIt = _inFileNames.begin(); flIt != _inFileNames.end(); ++flIt) {
		string locNam;
		for (size_t pos = 0; pos < flIt->size(); pos++) {
			if ( ((*flIt)[pos] == '/') || ((*flIt)[pos] == '\\') ) { // strip out *nix or Windows directory markers
				locNam.clear();
				continue;
			}
			if ( ((*flIt)[pos] == '_') || ((*flIt)[pos] == '.') ) {
				break;
			}
			locNam += (*flIt)[pos];
		}
		_lineNames.push_back(locNam);
	}
	
}

SFparse::SFparse(const string &fileList, const string &outFlNam, const unsigned long &alloc) : _outFileName(outFlNam), _bufAlloc(alloc) {
	string ext;
	bool foundDot = false;
	for (size_t pos = _outFileName.size() - 1; pos > 0; pos--) {
		if (foundDot) {
			if (_outFileName[pos] == '_') {
				break;
			}
			_chromName = _outFileName[pos] + _chromName;
			
		}
		if (_outFileName[pos] == '.') {
			foundDot = true;
		}
		if (!foundDot) {
			ext = _outFileName[pos] + ext;
		}
	}
	if (!foundDot) {
		cerr << "ERROR: no extension found in file " << _outFileName << " in SFparse extension-based constructor" << endl;
		exit(2);
	}
	if (_chromName.empty()) {
		cerr << "WARNING: no chromosome name found in output file name " << _outFileName << " in SFsparse constructor; setting default" << endl;
		_chromName = "NN";
		_chromNum  = 0;
	} else {
		_chromNum = 1;
	}
	if (ext == "bvt") {
		_outFileType = "BVT";
	} else if (ext == "bed") {
		_outFileType = "BED";
	} else {
		cerr << "ERROR: unknown extension " << ext << " for output file in SFparse extension-based constructor" << endl;
		exit(3);
	}
	ext.clear();
	auto outIt = _outFileName.end();
	_outFileName.erase(outIt - 4, outIt); // erase the extension
	
	ifstream lstIn(fileList.c_str());
	
	if (!lstIn) {
		cerr << "ERROR: cannot open file " << fileList << " listing files to process in SFsparse extension-based constructor" << endl;
		exit(1);
	}
	string eachFile;
	getline(lstIn, eachFile);
	foundDot = false;
	for (size_t pos = eachFile.size() - 1; pos > 0; pos--) {
		if (eachFile[pos] == '.') {
			foundDot = true;
			break;
		}
		ext = eachFile[pos] + ext;
	}
	if (ext == "seq") {
		_inFileType = "SEQ";
	} else {
		cerr << "ERROR: unknown extension " << ext << " for input files in SFparse extension-based constructor" << endl;
		exit(3);
	}
	
	if ( (eachFile[0] == 'r') && (eachFile[1] == ':') ) {
		eachFile.erase(0, 2);
		_refFlName = eachFile;
	} else {
		_inFileNames.push_back(eachFile);
	}
	while (getline(lstIn, eachFile)) {
		if ( (eachFile[0] == 'r') && (eachFile[1] == ':') ) {
			eachFile.erase(0, 2);
			_refFlName = eachFile;
		} else {
			_inFileNames.push_back(eachFile);
		}
	}
	
	lstIn.close();
	
	for (auto flIt = _inFileNames.begin(); flIt != _inFileNames.end(); ++flIt) {
		string locNam;
		for (size_t pos = 0; pos < flIt->size(); pos++) {
			if ( ((*flIt)[pos] == '/') || ((*flIt)[pos] == '\\') ) { // strip out *nix or Windows directory markers
				locNam.clear();
				continue;
			}
			if ( ((*flIt)[pos] == '_') || ((*flIt)[pos] == '.') ) {
				break;
			}
			locNam += (*flIt)[pos];
		}
		_lineNames.push_back(locNam);
	}
	
	if (_inFileNames.size() > numeric_limits<unsigned int>::max()) {
		cerr << "WARNING: number of lines " << _inFileNames.size() << " larger than allowed (" << numeric_limits<unsigned int>::max() << ")" << endl;
	}
}
SFparse& SFparse::operator=(const SFparse &inObj){
	if (this != &inObj) {
		_inFileNames = inObj._inFileNames;
		_lineNames   = inObj._lineNames;
		_refFlName   = inObj._refFlName;
		_outFileName = inObj._outFileName;
		_inFileType  = inObj._inFileType;
		_outFileType = inObj._outFileType;
		_chromName   = inObj._chromName;
		_chromNum    = inObj._chromNum;
		_bufAlloc    = inObj._bufAlloc;
		
	}
	
	return *this;
}

SFparse& SFparse::operator=(SFparse &&inObj){
	if (this != &inObj) {
		_inFileNames = move(inObj._inFileNames);
		_lineNames   = move(inObj._lineNames);
		_refFlName   = move(inObj._refFlName);
		_outFileName = move(inObj._outFileName);
		_inFileType  = move(inObj._inFileType);
		_outFileType = move(inObj._outFileType);
		_chromName   = move(inObj._chromName);
		_chromNum    = move(inObj._chromNum);
		_bufAlloc    = move(inObj._bufAlloc);
		
	}
	
	return *this;
}

void SFparse::operator()(){
	if ( (_inFileType == "SEQ") && (_outFileType == "BVT") ) {
		size_t bufSize = _bufAlloc/(_inFileNames.size() + 1);
		bool notDone = true;
		
		string fullOutName   = _outFileName + ".bvt";
		string outMetaFlName = _outFileName + ".bvtm";
		ofstream outMeta(outMetaFlName);
		if (!outMeta) {
			cerr << "ERROR: unable to open file " << outMetaFlName << " for metadata output in SFparse()" << endl;
			exit(6);
		}
		outMeta << _chromName << flush;
		for (auto lnNamIt = _lineNames.begin(); lnNamIt != _lineNames.end(); ++lnNamIt) {
			outMeta << " " << *lnNamIt << flush;
		}
		outMeta << endl;
		outMeta.close();
		
		ifstream inRef;
		
		vector<ifstream> inSeqs(_inFileNames.size());
		vector<char*> seqBufs(_inFileNames.size());
		
		remove(fullOutName.c_str());
		ofstream outDat(fullOutName, ios::binary);
		if (!outDat) {
			cerr << "ERROR: unable to open file " << fullOutName << " for data output in SFparse()" << endl;
			exit(6);
		}
		
		/*
		 *  There is a limit on how many files can be open at the same time
		 *  I have to close the SEQ files after reading each chunk; endPos is the place I save where I am to return to in the next iteration (if any)
		 */
		size_t endPosR      = 0;
		size_t endPosS      = 0;
		unsigned int chrPos = 1;
		
		// Read the FASTA files into the buffers, iterate until end of file is reached in the reference (this means that if, contrary to expectation, the sample files are longer they will be truncated)
		while (notDone) {
			char *refBuf = new char[bufSize]; // one extra for the null terminator
			inRef.open(_refFlName.c_str());
			if (!inRef) {
				cerr << "ERROR: unable to open reference file " << _refFlName << " in SFparse()" << endl;
				exit(5);
				
			}
			
			if (endPosR) {
				inRef.seekg(endPosR);
			}
			
			inRef.get(refBuf, bufSize);
			if (inRef.gcount() < bufSize - 1) { // did we read to the end?
				bufSize = inRef.gcount() + 1;
				notDone = false;
			}
			endPosS = endPosR;       // save the previous state of endPosR to read the population sample file
			endPosR = inRef.tellg(); // save position
			inRef.close();
			
			auto isfIt = inSeqs.begin();
			auto flnIt = _inFileNames.begin();
			for (auto sqbIt = seqBufs.begin(); sqbIt != seqBufs.end(); ++sqbIt) {
				isfIt->open(flnIt->c_str());
				if (!isfIt->is_open()) {
					cerr << "ERROR: unable to open file " << *flnIt << " in SFparse()" << endl;
					exit(5);
				}
				if (endPosS) {
					isfIt->seekg(endPosS);
				}
				
				*sqbIt = new char[bufSize];
				isfIt->get(*sqbIt, bufSize);
				isfIt->close();
				++isfIt;
				++flnIt;
			}
			
			char *polyLine = new char[_inFileNames.size() + 1];
			for (size_t i = 0; i < (bufSize - 1); i++) {
				bool polymorphic = false;
				polyLine[0] = refBuf[i];
				unsigned int iLine = 1;
				refBuf[i] = seqBufs[0][i]; // only looking for sites polymorphic within the sample; ones only divergent from reference not counted; therefore, the genotype of the i-th nucleotide for the first line is set to reference
				for (auto sbIt = seqBufs.begin(); sbIt != seqBufs.end(); sbIt++) {
					polyLine[iLine] = (*sbIt)[i];
					iLine++;
					if (refBuf[i] == 'N') {
						refBuf[i] = (*sbIt)[i];
					}
					if ( ((*sbIt)[i] != 'N') && ((*sbIt)[i] != refBuf[i]) ) { // if reference was missing, polymorphic definitely not set to true for this line
						polymorphic = true;
					}
				}
				if (polymorphic) {
					outDat.write(reinterpret_cast<char*>(&chrPos), sizeof(unsigned int));
					outDat.write(polyLine, _inFileNames.size() + 1);
				}
				chrPos++;
			}
			
			delete [] polyLine;
			
			for (auto sbIt = seqBufs.begin(); sbIt != seqBufs.end(); ++sbIt) {
				delete [] *sbIt;
			}
			delete [] refBuf;
		}
		
		outDat.close();
		
	} else if ( (_inFileType == "SEQ") && (_outFileType == "BED") ) {
		size_t bufSize = _bufAlloc/(_inFileNames.size() + 1);
		bool notDone = true;
		
		string outBedName = _outFileName + ".bed";
		string outBimName = _outFileName + ".bim";
		string outFamName = _outFileName + ".fam";
		
		// first save the .fam file
		ofstream outFam(outFamName);
		if (!outFam) {
			cerr << "ERROR: unable to open .fam file " << outFamName << " for output in SFparse()" << endl;
			exit(6);
		}
		for (auto lnNamIt = _lineNames.begin(); lnNamIt != _lineNames.end(); ++lnNamIt) {
			outFam << *lnNamIt << " " << *lnNamIt << " 0 0 0 -9" << endl;
		}
		outFam.close();
		
		ifstream inRef;
		
		vector<ifstream> inSeqs(_inFileNames.size());
		vector<char*> seqBufs(_inFileNames.size());
		
		remove(outBedName.c_str());
		ofstream outBed(outBedName, ios::binary);
		if (!outBed) {
			cerr << "ERROR: unable to open BED file " << outBedName << " for data output in SFparse()" << endl;
			exit(6);
		}
		
		remove(outBimName.c_str());
		ofstream outBim(outBimName);
		if (!outBim) {
			cerr << "ERROR: unable to open .bim file " << outBimName << " for data output in SFparse()" << endl;
			exit(6);
		}
		char magicBytes[] = {0x6C, 0x1B, 0x1}; // BED magic numbers go in the beginning of the file
		outBed.write(magicBytes, 3);
		
		/*
		 *  There is a limit on how many files can be open at the same time
		 *  I have to close the SEQ files after reading each chunk; endPos is the place I save where I am to return to in the next iteration (if any)
		 */
		size_t endPosR          = 0;
		size_t endPosS          = 0;
		unsigned int chrPos     = 1;
		unsigned int bedLineLen = ceil(static_cast<double>(_lineNames.size())/4.0); // SNPs are packed into bytes, four per byte with padding at each locus
		
		// Read the FASTA files into the buffers, iterate until end of file is reached in the reference (this means that if, contrary to expectation, the sample files are longer they will be truncated)
		while (notDone) {
			char *refBuf = new char[bufSize]; // one extra for the null terminator
			inRef.open(_refFlName.c_str());
			if (!inRef) {
				cerr << "ERROR: unable to open reference file " << _refFlName << " in SFparse()" << endl;
				exit(5);
				
			}
			
			if (endPosR) {
				inRef.seekg(endPosR);
			}
			
			inRef.get(refBuf, bufSize);
			if (inRef.gcount() < bufSize - 1) { // did we read to the end?
				bufSize = inRef.gcount() + 1;
				notDone = false;
			}
			endPosS = endPosR;       // save the previous state of endPosR to read the population sample file
			endPosR = inRef.tellg(); // save position
			inRef.close();
			
			auto isfIt = inSeqs.begin();
			auto flnIt = _inFileNames.begin();
			for (auto sqbIt = seqBufs.begin(); sqbIt != seqBufs.end(); ++sqbIt) {
				isfIt->open(flnIt->c_str());
				if (!isfIt->is_open()) {
					cerr << "ERROR: unable to open file " << *flnIt << " in SFparse()" << endl;
					exit(5);
				}
				if (endPosS) {
					isfIt->seekg(endPosS);
				}
				
				*sqbIt = new char[bufSize];
				isfIt->get(*sqbIt, bufSize);
				isfIt->close();
				++isfIt;
				++flnIt;
			}
			
			char *polyLine = new char[_inFileNames.size()];
			char *bedLine  = new char[bedLineLen];
			
			// Pre-form bit masks for going from a char array of genotypes to the packed BED format; the positions go in the reverse direction
			unordered_map<char, string> bitMasks;
			bitMasks['A'] = {static_cast<char>(0xFC), static_cast<char>(0xF3), static_cast<char>(0xCF), static_cast<char>(0x3F)}; // alternative (1/1 in plink)
			// no need for reference (2/2 in plink) because it will be all ones; doing it this way because SFS is heavy on low-frequency derived alleles and so I will mostly not have to do anything with bitmasks
			bitMasks['H'] = {static_cast<char>(0xFE), static_cast<char>(0xFB), static_cast<char>(0xEF), static_cast<char>(0xBF)}; // heterozygous; we actually do not have any so this is just for future development
			bitMasks['M'] = {static_cast<char>(0xFD), static_cast<char>(0xF7), static_cast<char>(0xDF), static_cast<char>(0x7F)}; // missing
			bitMasks['P'] = {static_cast<char>(0x3F), static_cast<char>(0x0F), static_cast<char>(0x03)};                          // padding
			
			// going over each site in the buffer, checking for polymorphism
			for (size_t i = 0; i < (bufSize - 1); i++) {
				bool polymorphic = false;
				bool biallelic   = true;    // only biallelic SNPs allowed in BED files
				char anc = refBuf[i];       // save the ancestral state
				char alt = '\0';
				unsigned int iLine = 0;
				refBuf[i] = seqBufs[0][i];  // only looking for sites polymorphic within the sample; ones only divergent from reference not counted; therefore, the genotype of the i-th nucleotide for the first line is set to reference
				
				// going over all the population lines
				for (auto sbIt = seqBufs.begin(); sbIt != seqBufs.end(); sbIt++) {
					polyLine[iLine] = (*sbIt)[i];
					iLine++;
					if (refBuf[i] == 'N') {
						refBuf[i] = (*sbIt)[i]; // this will keep happening until we hit a non-missing genotype
					}
					if ( ((*sbIt)[i] != 'N') && ((*sbIt)[i] != refBuf[i]) ) { // if reference was missing as of previous line, polymorphic definitely not set to true for this line because in that case we just set refBuf[i] to (*sbIt)[i] (that's why no else clause here!)
						if (!alt) {
							alt = (*sbIt)[i];
						} else {
							if (alt != (*sbIt)[i]) { // there already is an alternative and it is not the same as the current SNP
								biallelic = false;
								break; // if not biallelic, no use continuing with this site
							}
						}
						polymorphic = true;
					}
				}
				if (polymorphic && biallelic) {  // save a biallelic polymorphic site
					
					// first save the .bim metadata
					if (anc == 'N') { // label the SNP name with 'm' at the end is the ancestral state is missing
						outBim << _chromNum << " s" << chrPos << "m_" << _chromName << " -9 " << chrPos << " " << alt << " " << refBuf[i] << endl;
					} else if ( (alt != anc) && (refBuf[i] != anc) ) { // the SNP is biallelic in the sample, but the ancestral state is different from both
						outBim << _chromNum << " s" << chrPos << "d_" << _chromName << " -9 " << chrPos << " " << alt << " " << refBuf[i] << endl;
					} else {
						alt = (alt == anc ? refBuf[i] : alt); // assign ref to alt if alt is ancestral
						outBim << _chromNum << " s" << chrPos << "_" << _chromName << " -9 " << chrPos << " " << alt << " " << anc << endl;
					}
					
					unsigned int remainPad = bedLineLen * 4; // tracks how many genotypes are left in the padded line
					size_t iGeno = 0;                        // tracks the number of genotypes processed in the char array formed above (it is unpadded)
					for (size_t iBed = 0; iBed < bedLineLen; iBed++) {
						bedLine[iBed] = 0xFF;
						
						for (unsigned short bytePos = 0; bytePos < 4; bytePos++) { // actually going from the end of the byte, but that is already accounted for in the bitMaps object
							if (polyLine[iGeno] == alt) { // alternative (derived)
								bedLine[iBed] = bedLine[iBed] & bitMasks['A'][bytePos];
							} else if (polyLine[iGeno] == 'N') { // missing
								bedLine[iBed] = bedLine[iBed] & bitMasks['M'][bytePos];
							} // otherwise, it is reference and we do not do anything; no heterozygotes
							iGeno++;
							remainPad--;
							if (iGeno == _lineNames.size()) {
								if (remainPad != 0) {
									bedLine[iBed] = bedLine[iBed] & bitMasks['P'][remainPad - 1];
								}
								break; // that should automatically get us to the end of the outer loop, too
							}
							
						}
						
						
						
					}
					outBed.write(bedLine, bedLineLen);
				}
				chrPos++;
			}
			
			delete [] polyLine;
			delete [] bedLine;
			
			for (auto sbIt = seqBufs.begin(); sbIt != seqBufs.end(); ++sbIt) {
				delete [] *sbIt;
			}
			delete [] refBuf;
		}
		
		outBed.close();
		outBim.close();
		
		
		
	} else {
		cerr << "ERROR: unknown input or output format for parsing" << endl;
		exit(4);
	}
	
}


