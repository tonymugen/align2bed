
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
 * Class definitions and interface documentation for facilities that deal with common sequence and variant table file types.
 *
 */


#ifndef sequence_hpp
#define sequence_hpp

#include <vector>
#include <string>

using std::vector;
using std::string;
using std::move;

class SFparse;

/** \brief Sequence file parsing class
 *
 * Takes a list of files in one format and outputs one or more files in a different format, depending on settings. The data are presumed to come from a single chromosome.
 *
 * \note Formats will be added as the need arises.
 *
 */
class SFparse {
private:
	/// Input file names
	vector<string> _inFileNames;
	/// Line names
	vector<string> _lineNames;
	/// Reference file name
	string _refFlName;
	/// Output file name minus extention
	string _outFileName;
	/** \brief Input file format
	 *
	 * Supported input formats:
	 *
	 * - Headerless FASTA. Used in the DPGP project. Default extension is _.seq_
	 *
	 */
	string _inFileType;
	/** \brief Output file format
	 *
	 * Supported formats:
	 *
	 * - My own binary variant table. Default extension is _.bvt_ (Binary Variant Table). Variants are in rows, columns are chromosome position, reference nucleotide, string of base IDs (A,T,G,C,N) with no spaces. It is assumed that each chromosome is in a separate file. The file comes with a corresponding _bvtm_ file that has the metadata: chromosome name and names of the lines in a single space-separated line.
	 * - The _plink_ BED format. Default extension is _.bed_. It also comes with a _.bim_ and _.fam_ meta-data files.
	 *
	 */
	string _outFileType;
	/// Chromosome name
	string _chromName;
	/// Chromosome number
	unsigned short _chromNum;
	/** Buffer memory allocation 
	 *
	 * Total memory used for all input sequences. Default setting is 2000000000UL (2 Gb).
	 */
	unsigned long _bufAlloc;
	
public:
	/// Default constructor
	SFparse() : _bufAlloc(2000000000UL){};
	/** \brief Constructor with vectors of names
	 *
	 * Takes vectors of input and output file names. Note that the number of lines cannot be bigger than maximum of _unsigned int_. This is not checked. Also, the _lineNames_ vector must have one fewer elements than the _inFlNam_ vector.
	 *
	 * \param[in] inFlNam vector of input file names
	 * \param[in] lineNames vector of line names
	 * \param[in] refFlNam reference file name
	 * \param[in] outFlNam output file name
	 * \param[in] chrNam chromosome name
	 * \param[in] chrNum chromosome number
	 * \param[in] inFlType input file type
	 * \param[in] outFlType outout file type
	 * \param[in] alloc buffer allocation in bytes
	 */
	SFparse(const vector<string> &inFlNam, const vector<string> &lineNames, const string &refFlNam, const string &outFlNam, const string &chrNam, const unsigned short &chrNum, const string &inFlType, const string &outFlType, const unsigned long &alloc = 2000000000UL);
	/** \brief Constructor with file list and explicit types
	 *
	 * Takes a name of a file with the list of input files. The files should be listed one per line. Types to convert from and to are specified. The reference should be marked _r:_ in the list of files.
	 * Number of lines is checked and warning issued if it is larger than the _unsigned int_ maximum.
	 * Line names are derived from the non-reference files which should have the names before the underscore or extension (i.e., lineName_XXX.ext or lineName.ext). Directory names will be stripped out.
	 *
	 * \param[in] fileList name of the control file
	 * \param[in] outFlNam output file name
	 * \param[in] chrNam chromosome name
	 * \param[in] chrNum chromosome number
	 * \param[in] inFlType input file type
	 * \param[in] outFlType outout file type
	 * \param[in] alloc buffer allocation in bytes
	 */
	SFparse(const string &fileList, const string &outFlNam, const string &chrNam, const unsigned short &chrNum, const string &inFlType, const string &outFlType, const unsigned long &alloc = 2000000000UL);
	/** \brief Constructor with file list
	 *
	 * Takes a name of a file with the list of input files. The files should be listed one per line. Types are identified from file extensions. Only the first file listed for input is used (and checked) to determine the input type.
	 * The reference should be marked _r:_ in the list of files. Chromosome name is derived from the output file name which should immediately precede the extension and may be preceded by an underscore (i.e., XXX_chrName.ext).
	 * Number of lines is checked and warning issued if it is larger than the _unsigned int_ maximum. Line names are derived from the non-reference files which should have the names before the underscore or extension (i.e., lineName_XXX.ext or lineName.ext). Directory names will be stripped out.
	 *
	 * \param[in] fileList name of the control file
	 * \param[in] outFlNam output file name
	 * \param[in] alloc buffer allocation in bytes
	 */
	SFparse(const string &fileList, const string &outFlNam, const unsigned long &alloc = 2000000000UL);
	
	/// Destructor
	~SFparse(){};
	
	/** \brief Copy constructor
	 *
	 * \param[in] inObj object to be copied
	 */
	SFparse(const SFparse &inObj) : _inFileNames(inObj._inFileNames), _lineNames(inObj._lineNames), _refFlName(inObj._refFlName), _outFileName(inObj._outFileName), _inFileType(inObj._inFileType), _outFileType(inObj._outFileType), _chromName(inObj._chromName), _chromNum(inObj._chromNum), _bufAlloc(inObj._bufAlloc) {};
	/** \brief Copy assignement operator
	 *
	 * \param[in] inObj object to be copied
	 *
	 * \return SFparse object
	 */
	SFparse& operator=(const SFparse &inObj);
	/** \brief Move constructor
	 *
	 * \param[in] inObj object to be moved
	 */
	SFparse(SFparse &&inObj) : _inFileNames(move(inObj._inFileNames)), _lineNames(move(inObj._lineNames)), _refFlName(move(inObj._refFlName)), _outFileName(move(inObj._outFileName)), _inFileType(move(inObj._inFileType)), _outFileType(move(inObj._outFileType)), _chromName(move(inObj._chromName)), _chromNum(move(inObj._chromNum)), _bufAlloc(move(inObj._bufAlloc)) {};
	/** \brief Move assignement operator
	 *
	 * \param[in] inObj object to be moved
	 *
	 * \return SFparse object
	 */
	SFparse& operator=(SFparse &&inObj);
	
	/** \brief Modify output type
	 *
	 * \param[in] newType new output format
	 */
	void changeOutType(const string &newType) {_outFileType = newType; };
	
	/** \brief Input file parsing
	 *
	 * Function operator performs parsing of input files. Saves the results to the output file(s).
	 *
	 */
	void operator()();
	
};



#endif /* sequence_hpp */
