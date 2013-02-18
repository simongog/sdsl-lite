#ifndef CSA_TYPEDEFS
#define CSA_TYPEDEFS

#ifndef S_SA // if the suffix array sampling rate was not set
#define S_SA 64000 // set to default
#endif

#ifndef S_ISA // if the inverse suffix array sampling rate was not set
#define S_ISA 64000000 // set to default
#endif

#ifdef CSA_SADA
	#define CSA_TYPE csa_sada<enc_vector<coder::elias_delta,128>,S_SA,S_ISA>
	#define SUF "CSA_SADA"
#endif
#ifdef FM_HUFF_RRR255
	#define CSA_TYPE csa_wt<wt_huff<rrr_vector<255> >,S_SA,S_ISA>
	#define SUF "FM_HUFF_RRR255"
#endif
#ifdef FM_RLMN
	#define CSA_TYPE csa_wt<wt_rlmn<>,S_SA,S_ISA>
	#define SUF "FM_RLMN"
#endif
#ifdef FM_HUFF
	#define CSA_TYPE csa_wt<wt_huff<bit_vector,rank_support_v5<>,select_support_dummy,select_support_dummy>,S_SA,S_ISA>
	#define SUF "FM_HUFF"
#endif
#ifdef FM_HUFF_RRR127
	#define CSA_TYPE csa_wt<wt_huff<rrr_vector<127> >,S_SA,S_ISA>
	#define SUF "FM_HUFF_RRR127"
#endif
#ifdef FM_HUFF_RRR63
	#define CSA_TYPE csa_wt<wt_huff<rrr_vector<63> >,S_SA,S_ISA>
	#define SUF "FM_HUFF_RRR63"
#endif
#ifdef FM_HUFF_RRR15
	#define CSA_TYPE csa_wt<wt_huff<rrr_vector<15> >,S_SA,S_ISA>
	#define SUF "FM_HUFF_RRR15"
#endif

// set default CSA type
#ifndef CSA_TYPE
	#define INDEX CSA_SADA
	#define SUF "CSA_SADA"
#endif

#endif
