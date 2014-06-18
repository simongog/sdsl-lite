#include <iostream>
#include <limits>
#include <sdsl/vectors.hpp>
#include <sdsl/coder.hpp>
 
/**** Benchmark for self - delimiting codes ***********************************
For information about usage of this benchmark, see displayUsage - function.

To compile this benchmark, the following macros have to be defined 
(e.g. by passing them to compiler):
- VTYPES: a comma - separated list of sdsl vector types to be testet,
	e.g. vlc_vector<coder::elias_gamma>,vlc_vector<coder::elias_delta>
- VNAMES: symbolic names of the corresponding vector types, in same order 
	as in macro VTYPES, defined as a character array.
	According to the upper sample on macro VTYPES, VNAMES could be defined as
	{"VLC Vector with Elias Gamma Coder","VLC Vector with Elias Delta Coder"}
*/


//assert that needed macros are defined
#ifndef VTYPES
#error "Macro VTYPES with comma - separated list of vector types has to be \
	defined for compiling benchmark"
#endif

#ifndef VNAMES
#error "Macro VNAMES with an array of characters has to be \
	defined for compiling benchmark"
#endif

using namespace std;
using namespace sdsl;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

const char *(vectornames[]) = VNAMES;
const size_t vectorcount = sizeof(vectornames) / sizeof(vectornames[0]);

struct iv_testresult { //testcase for one defined int vector
	double enc_MBperSec; //encoding rate: megabytes per second
	double dec_MBperSec; //decoding rate: megabytes per second
	double comp_percent; //compression rate: needed space in percentage compared
			//to original integer vector
};

//benchmark method declaration
template<class... Vectors> //used vectors for benchmark
bool runTestcase( const int_vector<> &iv, iv_testresult *result );

//stuff for nice printing
void displayUsage(const char *pname);
void displayHeading();
void displayResult( const char *testcase, const iv_testresult *result );

int main(const int argc, const char **argv)
{
	//check args
	if ((argc - 1) % 3 != 0) {
		displayUsage(argv[0]);
		return 1;
	}

	//set up needed structures
	const size_t testcasecount = (argc - 1) / 3;
	iv_testresult overallresult[vectorcount];
	
	//prepare overall result
	for (size_t i = 0; i < vectorcount; i++) {
		overallresult[i].enc_MBperSec = 0.0;
		overallresult[i].dec_MBperSec = 0.0;
		overallresult[i].comp_percent = 0.0;
	}

	//start fetching test cases and run benchmark
	displayHeading();
	for (size_t i = 0; i < testcasecount; i++) {
		const char *testcase = argv[3*i + 1];
		const char *file =     argv[3*i + 2]; //file of saved vector
		const char *type =     argv[3*i + 3]; //type of saved vector
		uint8_t v_type = type[0]=='d' ? 'd' : type[0] - '0';

		//load vector
		int_vector<> iv;
		if (!load_vector_from_file(iv, file, v_type)) {
			cerr << "ERROR: vector from file " << file 
			     << " with type " << type << " could not be loaded" 
			     << endl;
			displayUsage(argv[0]);
			return 1;
		}

		//run test
		iv_testresult result[vectorcount];
		if (!runTestcase<VTYPES>( iv, result )) {
			cerr << "Testcase " << testcase << "failed" << endl;
			return 1;
		}

		//print result
		displayResult( testcase, result );

		//and sum up results for overall result
		for (size_t j = 0; j < vectorcount; j++) {
			overallresult[j].enc_MBperSec += result[j].enc_MBperSec;
			overallresult[j].dec_MBperSec += result[j].dec_MBperSec;
			overallresult[j].comp_percent += result[j].comp_percent;
		}
	}

	//build average for overall result
	for (size_t i = 0; i < vectorcount; i++) {
		overallresult[i].enc_MBperSec /= testcasecount;
		overallresult[i].dec_MBperSec /= testcasecount;
		overallresult[i].comp_percent /= testcasecount;
	}

	//and display overall results
	displayResult( "Overall", overallresult );
	return 0;
}

//// BENCHMARK METHODS ////////////////////////////////////////////////////////
template<class Vector> //used compression vector type
bool runSingleTest( const int_vector<> &testcase, iv_testresult &result ) {
	//test encoding rate by constructing Vector
	auto start = timer::now();
	Vector test( testcase );
	auto stop = timer::now();
	result.enc_MBperSec = size_in_mega_bytes( testcase )
		 / duration_cast<seconds>(stop-start).count();

	//care for compression rate
	result.comp_percent = size_in_mega_bytes(test) 
		/ size_in_mega_bytes(testcase) * 100.0;

	//and finally for decoding rate
	//use a trick to decode all values: since (currently) all vectors are
	//using sample tables, access the element right before the next sampling
	//entry, so everything between 2 samples has to be decoded.
	size_t sample_dens = test.get_sample_dens();
	start = timer::now();
	//repeat test 5 times to avoid infinite decoding rates
	for (size_t j = 0; j < 5; j++) {
		size_t i = sample_dens - 1;
		for (; i < test.size(); i += sample_dens) {
			test[i]; //acess element right before next sample entry
		}
		//and finally access last element if not done yet
		if (i != test.size() + sample_dens - 1)
			test[test.size() - 1];
	}
	stop = timer::now();
	result.dec_MBperSec = size_in_mega_bytes( testcase )
		 / duration_cast<seconds>(stop-start).count() 
		* 5.0; //multiply with 5 since vector was decoded 5 times

	return true; //may use this return type for error detection in future
}

template<class... Vectors> //used vectors for benchmark
bool runTestcase( const int_vector<> &testcase, iv_testresult *result ) {
	size_t i = 0;
	//do variadic template pack expansion
	bool testfine[] = { runSingleTest<Vectors>( testcase, result[i++] )... };
	bool testsfine = true;
	for (i = 0; i < vectorcount; i++) {
		if (!testfine[i]) {
			cerr << "Test on Vector " << vectornames[i]
			     << "failed" << endl;
			testsfine = false;
		}
	}
	return testsfine;	
}

//// DISPLAYING OF RESULTS ////////////////////////////////////////////////////

void displayUsage(const char *pname) {
	cerr << "USAGE: " << pname << " [testcase file vectortype]*"
	     << endl;
	cerr << "DESCRIPTION:" << endl;
	cerr << "\tThis Program runs a benchmark on self-delimiting "
                    << "Codes." << endl;
	cerr << "\tProgram needs triples of parameters "
	     << "for each test case, see Parameter section." << endl;
	cerr << "\tProgram will test a couple of compression vectors " 
	     << endl
	     << "\ton measured encoding and decoding rates," << endl
	     << "\tplus the compression rate in percent " 
	     << "(compared to the original integer vector)" << endl
	     << "\tfor each testcase." 
	     << endl
	     << "\tAdditionally, an overall result on different " 
	     << endl << "\tcompression vectors is printed." << endl;
	cerr << "\tThe generated output uses a CSV format, so " 
	     << "you may save it to a csv file for better visability"
	     << endl << "\tand other utilites." << endl;
	cerr << "PARAMETERS: The parameters have to be passed as "
	     << " triples for each test case." << endl
	     << "\tA Triple consist of " << endl
	     << "\t\t- testcase: A name for the test case" << endl
	     << "\t\t- file: a path to the file where the test case" << endl
	     << "\t\t\t(an integer vector) is contained" << endl
	     << "\t\t- vectortype: type of saved integer vector" << endl
	     << "\t\t\t0: serialized int_vector<>" << endl
	     << "\t\t\t1: byte sequence" << endl
	     << "\t\t\t2: 16-bit word sequence" << endl
	     << "\t\t\t4: 32-bit word sequence" << endl
	     << "\t\t\t8: 64-bit word sequence" << endl
	     << "\t\t\td: Parse decimal numbers" << endl;
	cerr << "TESTET COMPRESSION VECTORS:" << endl;
	for (size_t i = 0; i < vectorcount; i++) {
		cerr << "\t- " << vectornames[i] << endl;
	}
}
void displayHeading() {
	cout << left; //left justify
	//add a comment how to read values
	cout << "# encoding / decoding rate unit: MB/s" << endl;
	cout << "# compression : percentage of needed space "
	     << " compared to original vector" << endl;
	//and print a header for csv output
	cout << setw(20) << "testcase" 
	     << setw(1) << ";" << setw(20) << "vector"
	     << setw(1) << ";" << setw(20) << "encodingrate"
	     << setw(1) << ";" << setw(20) << "decodingrate"
	     << setw(1) << ";" << "compressionrate" << endl;
}

void displayResult( const char *testcase, const iv_testresult *result ) {
	cout << left << fixed; //prepare cout
	for (size_t i = 0; i < vectorcount; i++) {
		cout << setw(20) << testcase 
		     << setw(1) << ";" << setw(20) << vectornames[i]
		     << setw(1) << ";" << setw(20) << result[i].enc_MBperSec
		     << setw(1) << ";" << setw(20) << result[i].dec_MBperSec
		     << setw(1) << ";" << result[i].comp_percent << endl;
	}
}
