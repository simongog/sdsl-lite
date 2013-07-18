/* sdsl - succinct data structures library
    Copyright (C) 2012-2013 Simon Gog
    Copyright (C) 2013 Timo Beller
*/

#include "sdsl/construct_lcp.hpp"
#include "sdsl/construct.hpp"
#include <stdexcept>
#include <algorithm>

namespace sdsl
{

void construct_lcp_semi_extern_PHI(cache_config& config)
{
    typedef int_vector<>::size_type size_type;
    int_vector_file_buffer<> sa_buf(cache_file_name(constants::KEY_SA, config));
    size_type n = sa_buf.int_vector_size;
    if (1==n) {
        int_vector<> lcp(1, 0);
        store_to_cache(lcp, constants::KEY_LCP, config);
        return;
    }
    const uint8_t log_q = 6; // => q=64
    const uint32_t q = 1<<log_q;
    const uint64_t modq = bits::lo_set[log_q];

    // n-1 is the maximum entry in SA
    int_vector<64> plcp((n-1+q)>>log_q);

    mm::log("lcp-calc-sparse-phi-begin");
    sa_buf.reset();
    for (size_type i=0, sai_1=0, r=0, r_sum=0; r_sum < n;) {   // we can start at i=0. if SA[i]%q==0
        for (; i < r_sum+r; ++i) {               // we set PHI[(SA[i]=n-1)%q]=0, since T[0]!=T[n-1]
            size_type sai = sa_buf[i-r_sum];
            if ((sai & modq) == 0) {
                if ((sai>>log_q) >= plcp.size()) {
//                    std::cerr<<"sai="<<sai<<" log_q="<<log_q<<" sai>>log_q="<<(sai>>log_q)<<" "<<sai_1<<std::endl;
//                    std::cerr<<"n="<<n<<" "<<" plcp.size()="<<plcp.size();
                }
                plcp[sai>>log_q] = sai_1;
            }
            sai_1 = sai;
        }
        r_sum += r; r = sa_buf.load_next_block();
    }
    mm::log("lcp-calc-sparse-phi-end");

    mm::log("lcp-load-text-begin");
    int_vector<8> text;
    load_from_cache(text, constants::KEY_TEXT, config);
    mm::log("lcp-load-text-end");


    mm::log("lcp-calc-sparse-plcp-begin");
    for (size_type i=0,j,k,l=0; i < plcp.size(); ++i) {
        j =	i<<log_q;   // j=i*q
        k = plcp[i];
        while (text[j+l] == text[k+l])
            ++l;
        plcp[i] = l;
        if (l >= q) {
            l -= q;
        } else {
            l = 0;
        }
    }
    mm::log("lcp-calc-sparse-plcp-end");

    osfstream lcp_out_buf(cache_file_name(constants::KEY_LCP, config), std::ios::binary | std::ios::trunc);   // open buffer for plcp
    size_type bit_size = n*sa_buf.width;
    lcp_out_buf.write((char*) &(bit_size), sizeof(sa_buf.int_vector_size));	// write size of vector
    lcp_out_buf.write((char*) &(sa_buf.width), sizeof(sa_buf.width));       // write int_width of vector
    size_type wb = 0;

    size_type buffer_size = 4000000; // buffer_size is a multiple of 8!
    sa_buf.reset(buffer_size);
    for (size_type i=0, r=0, r_sum=0, sai_1=0,l=0, sai=0,iq=0; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            /*size_type*/ sai = sa_buf[i-r_sum];
//				std::cerr<<"i="<<i<<" sai="<<sai<<std::endl;
            if ((sai & modq) == 0) { // we have already worked the value out ;)
                sa_buf.set_int(i-r_sum, l=plcp[sai>>log_q]);
            } else {
                /*size_type*/ iq = sai & bits::lo_unset[log_q];
                l  = plcp[sai>>log_q];
                if (l >	(sai-iq))
                    l -= (sai-iq);
                else
                    l=0;
                while (text[ sai+l ] == text[ sai_1+l ])
                    ++l;
                sa_buf.set_int(i-r_sum, l);
            }
#ifdef CHECK_LCP
            size_type j=0;
            for (j=0; j<l; ++j) {
                if (text[sai+j] !=text[sai_1+j]) {
                    std::cout<<"lcp["<<i<<"]="<<l<<" is two big! "<<j<<" is right!"<<" sai="<<sai<<std::endl;
                    if ((sai&modq)!=0)
                        std::cout<<" plcp[sai>>log_q]="<<plcp[sai>>log_q]<<" sai-iq="<<sai-iq<<" sai="<<sai<<" sai-iq="<<sai-iq<<std::endl;
                    break;
                }
            }
#endif
            sai_1 = sai;
        }
        if (r > 0) {
            size_type cur_wb = (r*sa_buf.width+7)/8;
            lcp_out_buf.write((const char*)sa_buf.data(), cur_wb);
            wb += cur_wb;
        }
        r_sum += r; r = sa_buf.load_next_block();
    }

    if (wb%8) {
        lcp_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
    }
    lcp_out_buf.close();
    register_cache_file(constants::KEY_LCP, config);
    return;
}

void construct_lcp_go(cache_config& config)
{
    typedef int_vector<>::size_type size_type;
#ifdef STUDY_INFORMATIONS
    size_type racs  = 0; // random accesses to the text
    size_type matches = 0;
    size_type comps2 = 0; // comparisons the second phase
#endif
    int_vector<8> text;
    load_from_cache(text, constants::KEY_TEXT, config);
    int_vector_file_buffer<> sa_buf(cache_file_name(constants::KEY_SA, config));   // initialize buffer for suffix array
    const size_type n = sa_buf.int_vector_size;
    const size_type m = 254; // LCP[i] == m+1 corresp. to LCP[i]>= m+1; LCP[i] <= m corresp. to LCP[i] was calculated

    if (1==n) {
        int_vector<> lcp(1, 0);
        store_to_cache(lcp, constants::KEY_LCP, config);
        return;
    }

    size_type cnt_c[257] = {0};   // counter for each character in the text
    size_type cnt_cc[257] = {0};  // prefix sum of the counter cnt_c
    size_type cnt_cc2[257] = {0};  //
    size_type omitted_c[257] = {0};  // counts the omitted occurrences for the second phase
    size_type prev_occ_in_bwt[256] = {0};  // position of the previous occurrence of each character c in the bwt
    for (size_type i=0; i<256; ++i) prev_occ_in_bwt[i] = (size_type)-1; // initialize the array with -1
    unsigned char alphabet[257] = {0};
    uint8_t sigma = 0;

    tLI m_list[2][256];
    size_type m_char_count[2] = {0};
    uint8_t m_chars[2][256] = {{0},{0}};

    size_type nn = 0; // n' for phase 2
    // phase 1: calculate lcp_sml;  memory consumption: 2n bytes (lcp_sml=n bytes, text=n bytes)
    {

        int_vector<8> lcp_sml(n, 0); // initialize array for small values of first phase; note lcp[0]=0
        size_type done_cnt=0;

        for (size_type i=0; i<n; ++i) { // initialize cnt_c
            ++cnt_c[text[i]+1];
        }
        for (int i=1; i<257; ++i) { // calculate sigma and initailize cnt_cc
            if (cnt_c[i] > 0) {
                alphabet[sigma++] = (unsigned char)(i-1);
            }
            cnt_cc[i] = cnt_c[i] + cnt_cc[i-1];
        }
        alphabet[sigma] = '\0';
        {
            int_vector_file_buffer<8> bwt_buf(cache_file_name(constants::KEY_BWT, config)); // initialize buffer of bwt
            sa_buf.reset();
            sa_buf.load_next_block();
            size_type sai_1 = sa_buf[0];  // store value of sa[i-1]
            size_type r = bwt_buf.load_next_block();
            uint8_t bwti_1 = bwt_buf[0];       // store value of BWT[i-1]
            lcp_sml[ cnt_cc[bwti_1]++ ] = 0;   // lcp_sml[ LF[0] ] = 0
            prev_occ_in_bwt[bwti_1] = 0;  // init previous occurence of character BWT[0]
            ++omitted_c[alphabet[0]];	  //

            int_vector<64> rmq_stack(2*(m+10));   // initialize stack for m+10 elements representing (position, value)
            rmq_stack[0] = 0; rmq_stack[1] = 0;  // first element (-1, -1)
            rmq_stack[2] = 1; rmq_stack[3] = 0;  // second element (0, -1)
            size_type rmq_end=3;				 // index of the value of the topmost element

            const size_type m_mod2 = m%2;
            uint8_t cur_c = alphabet[1];
            size_type big_val = 0;
            for (size_type i=1, sai, r_sum=0, cur_c_idx=1, cur_c_cnt=cnt_c[alphabet[1]+1]; r_sum < n;) {
                for (; i < r_sum+r; ++i, --cur_c_cnt) {
                    uint8_t bwti = bwt_buf[i-r_sum];
                    sai = sa_buf[i-r_sum];
                    size_type lf = cnt_cc[bwti];
                    if (!cur_c_cnt) {// cur_c_cnt==0, if there is no more occurence of the current character
                        if (cur_c_cnt < sigma) {
                            cur_c_cnt = cnt_c[(cur_c=alphabet[++cur_c_idx])+1];
                        }
                    }
                    size_type l=0;
                    if (i >= cnt_cc[cur_c]) { // if the current lcp entry is not already done   TODO: schleife von i bis cnt_cc[cur_c]
                        if (bwti == bwti_1 and lf < i) {  // BWT[i]==BWT[i-1]
                            l = lcp_sml[lf] ? lcp_sml[lf]-1 : 0; // l = LCP[LF[i]]-1; l < m+1
                            if (l == m) { // if LCP[LF[i]] == m+1; otherwise LCP[LF[i]] < m+1  the result is correct
                                l += (text[sai_1+m] == text[sai+m]);
#ifdef STUDY_INFORMATIONS
                                if ((sai_1^sai)>>6) // if i and phii are in the same cache line
                                    ++racs;
#endif
                            }
                            lcp_sml[i] = l;
                            ++done_cnt;
                        } else { // BWT[i] != BWT[i-1] or LF[i] > i
                            if (lf < i)
                                l = lcp_sml[lf] ? lcp_sml[lf]-1 : 0;
#ifdef STUDY_INFORMATIONS
                            if ((sai_1^sai)>>6) // if i and phii are in the same cache line
                                ++racs;
#endif
                            while (text[sai_1+l] == text[sai+l] and l < m+1) {
                                ++l;
#ifdef STUDY_INFORMATIONS
                                ++matches;
#endif
                            }
                            lcp_sml[i] = l;
                        }
                    } else { // if already done
                        l = lcp_sml[i];  // load LCP value
                    }
                    if (l > m) {
                        ++big_val;
                        if (i > 10000 and i < 10500 and big_val > 3000) { // if most of the values are big: switch to PHI algorithm
                            util::clear(text);
                            util::clear(lcp_sml);
                            construct_lcp_PHI<8>(config);
                            return;
                        }
                    }
                    // invariant: l <= m+1
                    // begin update rmq_stack
                    size_type x = l+1;
                    size_type j = rmq_end;
                    while (x <= rmq_stack[j]) j-=2;  // pop all elements with value >= l
                    rmq_stack[++j] = i+1; // push position i
                    rmq_stack[++j] = x;	  // push value	l
                    rmq_end = j;          // update index of the value of the topmost element
                    if (lf > i) {   // if LF[i] > i, we can calculate LCP[LF[i]] in constant time with rmq
                        ++done_cnt;
                        // rmq query for lcp-values in the interval I=[prev_occ_in_bwt[BWT[i]]+1..i]
                        // rmq is linear in the stack size; can also be implemented with binary search on the stack
                        size_type x_pos = prev_occ_in_bwt[bwti]+2;
                        size_type j = rmq_end-3;
                        while (x_pos <= rmq_stack[j]) j-=2;   //  search smallest value in the interval I
                        lcp_sml[lf] = rmq_stack[j+3] - (rmq_stack[j+3]==m+2); // if lcp-value equals m+1, we subtract 1
                    }
                    if (l >= m) {
                        if (l == m)
                            push_front_m_index(nn, cur_c, m_list[m_mod2], m_chars[m_mod2], m_char_count[m_mod2]);
                        ++nn;
                    } else
                        ++omitted_c[cur_c];

                    prev_occ_in_bwt[bwti] = i;  	 // update previous position information for character BWT[i]
                    ++cnt_cc[bwti];					 // update counter and therefore the LF information
                    sai_1 = sai;					 // update SA[i-1]
                    bwti_1 = bwti;					 // update BWT[i-1]
                }
                r_sum += r; r = bwt_buf.load_next_block();
                sa_buf.load_next_block();
            }
        }
        util::clear(text);

        if (n > 1000 and nn > 5*(n/6)) {  // if we would occupy more space than the PHI algorithm => switch to PHI algorithm
            util::clear(lcp_sml);
            construct_lcp_PHI<8>(config);
            return;
        }
        store_to_cache(lcp_sml, "lcp_sml", config);
    }
#ifdef STUDY_INFORMATIONS
    std::cout<<"# n="<<n<<" nn="<<nn<<" nn/n="<<((double)nn)/n<<std::endl;
#endif

    // phase 2: calculate lcp_big
    {
//			std::cout<<"# begin calculating LF' values"<<std::endl;
        int_vector<> lcp_big(nn, 0, bits::hi(n-1)+1); // lcp_big first contains adapted LF values and finally the big LCP values
        {
            // initialize lcp_big with adapted LF values
            bit_vector todo(n,0);  // bit_vector todo indicates which values are >= m in lcp_sml
            {
                // initialize bit_vector todo
                int_vector_file_buffer<8> lcp_sml_buf(cache_file_name("lcp_sml", config)); // load lcp_sml
                for (size_type i=0, r_sum=0, r = lcp_sml_buf.load_next_block(); r_sum < n;) {
                    for (; i < r_sum+r; ++i) {
                        if (lcp_sml_buf[i-r_sum] >= m) {
                            todo[i] = 1;
                        }
                    }
                    r_sum += r; r = lcp_sml_buf.load_next_block();
                }
            }

            cnt_cc2[0] = cnt_cc[0]= 0;
            for (size_type i=1, omitted_sum=0; i<257; ++i) { // initialize cnt_cc and cnt_cc2
                cnt_cc[i] = cnt_c[i] + cnt_cc[i-1];
                omitted_sum += omitted_c[i-1];
                cnt_cc2[i] = cnt_cc[i] - omitted_sum;
            }

            int_vector_file_buffer<8> bwt_buf(cache_file_name(constants::KEY_BWT, config)); // load BWT
            for (size_type i=0, i2=0, r_sum=0, r=bwt_buf.load_next_block(); r_sum < n;) {
                for (; i < r_sum+r; ++i) {
                    uint8_t b = bwt_buf[i - r_sum];  // store BWT[i]
                    size_type lf_i = cnt_cc[b]; // LF[i]
                    if (todo[i]) { // LCP[i] is a big value
                        if (todo[lf_i]) { // LCP[LF[i]] is a big entry
                            lcp_big[i2] = cnt_cc2[b]; // LF'[i]
                        }/*else{
								lcp_big[i2] = 0;
							}*/
                        ++i2;
                    }
                    if (todo[lf_i]) { // LCP[LF[i]] is a big entry
                        ++cnt_cc2[b];	// increment counter for adapted LF
                    }
                    ++cnt_cc[b]; // increment counter for LF
                }
                r_sum += r; r = bwt_buf.load_next_block();
            }
        }

//			std::cout<<"# begin initializing bwt2, shift_bwt2, run2"<<std::endl;
        int_vector<8> bwt2(nn), shift_bwt2(nn); // BWT of big LCP values, and shifted BWT of big LCP values
        bit_vector run2(nn+1);					// indicates for each entry i, if i and i-1 are both big LCP values
        run2[nn] = 0;							// index nn is not a big LCP value
        {
            // initialize bwt2, shift_bwt2, adj2
            int_vector_file_buffer<8> lcp_sml_buf(cache_file_name("lcp_sml", config)); // load lcp_sml
            int_vector_file_buffer<8> bwt_buf(cache_file_name(constants::KEY_BWT, config)); // load BWT
            uint8_t b_1 = '\0'; // BWT[i-1]
            bool is_run = false;
            for (size_type i=0, i2=0, r_sum=0, r = 0; r_sum < n;) {
                for (; i < r_sum+r; ++i) {
                    uint8_t b = bwt_buf[i-r_sum];
                    if (lcp_sml_buf[i-r_sum] >= m) {
                        bwt2[i2] 		= b;
                        shift_bwt2[i2]	= b_1;
                        run2[i2] 		= is_run;
                        is_run 			= true;
                        ++i2;
                    } else {
                        is_run = false;
                    }
                    b_1 = b;
                }
                r_sum += r; r = lcp_sml_buf.load_next_block();
                bwt_buf.load_next_block();
            }
        }

        bit_vector todo2(nn+1, 1); // init all values with 1, except
        todo2[nn] = 0; 			   // the last one! (handels case "i < nn")

//			std::cout<<"# begin calculating m-indices"<<std::endl;
        {
            // calculate m-indices, (m+1)-indices,... until we are done
            size_type m2 = m;
            size_type char_ex[256]; for (size_type i=0; i<256; ++i) char_ex[i] = nn;
            size_type char_occ=0;
            size_type m_mod2 = m2%2, mm1_mod2 = (m2+1)%2;
            while (m_char_count[m_mod2] > 0) { // while there are m-indices, calculate (m+1)-indices and write m-indices
                // For all values LCP[i] >= m2 it follows that todo2[i] == 1
                // list m_list[mm1_mod2][b] is sorted in decreasing order
                ++m2;
                mm1_mod2 = (m2+1)%2, m_mod2 = m2%2;
                m_char_count[m_mod2] = 0;

                std::sort(m_chars[mm1_mod2], m_chars[mm1_mod2]+m_char_count[mm1_mod2]); // TODO: ersetzen?

                for (size_type mc=0; mc<m_char_count[mm1_mod2]; ++mc) { // for every character
                    tLI& mm1_mc_list = m_list[mm1_mod2][m_chars[mm1_mod2][ m_char_count[mm1_mod2]-1-  mc ]];
//						size_type old_i = nn;
                    while (!mm1_mc_list.empty()) {
                        size_type i = mm1_mc_list.front();  // i in [0..n-1]
                        mm1_mc_list.pop_front();
                        // For all values LCP[i] >= m-1 it follows that todo2[i] == 1
                        for (size_type k=i; todo2[k]; --k) {
#ifdef STUDY_INFORMATIONS
                            ++comps2;
#endif
                            uint8_t b = shift_bwt2[k];
                            if (char_ex[b] != i) {
                                char_ex[b] = i;
                                ++char_occ;
                            }
                            if (!run2[k])
                                break;
                        }
                        for (size_type k=i; todo2[k] and char_occ; ++k) {
#ifdef STUDY_INFORMATIONS
                            ++comps2;
#endif
                            uint8_t b = bwt2[k];
                            if (char_ex[b] == i) {
                                size_type p = lcp_big[k];
                                push_back_m_index(p, b, m_list[m_mod2], m_chars[m_mod2], m_char_count[m_mod2]);
                                char_ex[b] = nn;
                                --char_occ;
                            }
                            if (!run2[k+1])
                                break;
                        }
                        lcp_big[ i ] = m2-1;
                        todo2[ i ] = 0;
//							old_i = i;
                    }
                }
            }

        }

        store_to_cache(lcp_big, "lcp_big", config);
    } // end phase 2

//		std::cout<<"# merge lcp_sml and lcp_big"<<std::endl;
    // phase 3: merge lcp_sml and lcp_big and save to disk
    {
        const size_type buffer_size = 1000000; // buffer_size has to be a multiple of 8!
        int_vector_file_buffer<> lcp_big_buf(cache_file_name("lcp_big", config)); 				// file buffer containing the big LCP values
        int_vector_file_buffer<8> lcp_sml_buf(cache_file_name("lcp_sml", config), buffer_size);// file buffer containing the small LCP values

        osfstream lcp_out_buf(cache_file_name(constants::KEY_LCP, config), std::ios::binary | std::ios::trunc);    // open out file stream
        int_vector<> lcp_buf(buffer_size, 0, lcp_big_buf.width); // buffer for the resulting LCP array
        size_type bit_size = n*lcp_big_buf.width;
        lcp_out_buf.write((char*) &(bit_size), sizeof(lcp_big_buf.int_vector_size));		   // write size
        lcp_out_buf.write((char*) &(lcp_big_buf.width), sizeof(lcp_big_buf.width));   // write int_width
        size_type wb = 0; // written bytes of data

        for (size_type i=0, r_sum=0, r=0, i2=0 ,r2_sum=0, r2=0; r_sum < n;) {
            for (; i < r_sum + r; ++i) {
                size_type l = lcp_sml_buf[i-r_sum];
                if (l >= m) { // if l >= m it is stored in lcp_big
                    if (i2 >= r2_sum + r2) {
                        r2_sum += r2; r2 = lcp_big_buf.load_next_block();
                    }
                    l = lcp_big_buf[i2-r2_sum];
                    ++i2;
                }
                lcp_buf[i-r_sum] = l;
            }
            if (r > 0) {
                size_type cur_wb = r*lcp_buf.width();
                if (cur_wb%8 != 0)
                    cur_wb += 7;
                cur_wb /= 8;
                lcp_out_buf.write((const char*)lcp_buf.data(), cur_wb);
                wb += cur_wb;
            }
            r_sum += r; r = lcp_sml_buf.load_next_block();
        }
        if (wb%8) {
            lcp_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
        }
    }
    register_cache_file(constants::KEY_LCP, config);

#ifdef STUDY_INFORMATIONS
    std::cout<<"# racs: "<<racs<<std::endl;
    std::cout<<"# matches: "<<matches<<std::endl;
    std::cout<<"# comps2: "<<comps2<<std::endl;
#endif
    return;
}

void construct_lcp_goPHI(cache_config& config)
{
    typedef int_vector<>::size_type size_type;
    int_vector<8> text;
    load_from_cache(text, constants::KEY_TEXT, config);  // load text from file system
    int_vector_file_buffer<> sa_buf(cache_file_name(constants::KEY_SA, config));   // initialize buffer for suffix array
    const size_type n = sa_buf.int_vector_size;
    const size_type m = 254; // LCP[i] == m+1 corresp. to LCP[i]>= m+1; LCP[i] <= m corresp. to LCP[i] was calculated

    if (1==n) {
        int_vector<> lcp(1, 0);
        store_to_cache(lcp, constants::KEY_LCP, config);
        return;
    }

    size_type cnt_c[257] = {0};   // counter for each character in the text
    size_type cnt_cc[257] = {0};  // prefix sum of the counter cnt_c
    size_type omitted_c[257] = {0};  // counts the omitted occurrences for the second phase
    size_type prev_occ_in_bwt[256] = {0};  // position of the previous occurrence of each character c in the bwt
    for (size_type i=0; i<256; ++i) prev_occ_in_bwt[i] = (size_type)-1; // initialize the array with -1
    unsigned char alphabet[257] = {0};
    uint8_t sigma = 0;

    size_type nn = 0; // n' for phase 2
    // phase 1: calculate lcp_sml;  memory consumption: 2n bytes (lcp_sml=n bytes, text=n bytes)
    {

        int_vector<8> lcp_sml(n, 0); // initialize array for small values of first phase; note lcp[0]=0
        size_type done_cnt=0;

        for (size_type i=0; i<n; ++i) { // initialize cnt_c
            ++cnt_c[text[i]+1];
        }
        for (int i=1; i<257; ++i) { // calculate sigma and initailize cnt_cc
            if (cnt_c[i] > 0) {
                alphabet[sigma++] = (unsigned char)(i-1);
            }
            cnt_cc[i] = cnt_c[i] + cnt_cc[i-1];
        }
        alphabet[sigma] = '\0';
        {
            int_vector_file_buffer<8> bwt_buf(cache_file_name(constants::KEY_BWT, config)); // initialize buffer of bwt
            sa_buf.reset();
            sa_buf.load_next_block();
            size_type sai_1 = sa_buf[0];  // store value of sa[i-1]
            size_type r = bwt_buf.load_next_block();
            uint8_t bwti_1 = bwt_buf[0];       // store value of BWT[i-1]
            lcp_sml[ cnt_cc[bwti_1]++ ] = 0;   // lcp_sml[ LF[0] ] = 0
            prev_occ_in_bwt[bwti_1] = 0;  // init previous occurence of character BWT[0]
            ++omitted_c[alphabet[0]];	  //

            int_vector<64> rmq_stack(2*(m+10));   // initialize stack for m+10 elements representing (position, value)
            rmq_stack[0] = 0; rmq_stack[1] = 0;  // first element (-1, -1)
            rmq_stack[2] = 1; rmq_stack[3] = 0;  // second element (0, -1)
            size_type rmq_end=3;				 // index of the value of the topmost element

            uint8_t cur_c = alphabet[1];
            for (size_type i=1, sai, r_sum=0, cur_c_idx=1, cur_c_cnt=cnt_c[alphabet[1]+1]; r_sum < n;) {
                for (; i < r_sum+r; ++i, --cur_c_cnt) {
                    uint8_t bwti = bwt_buf[i-r_sum];
                    sai = sa_buf[i-r_sum];
                    size_type lf = cnt_cc[bwti];
                    if (!cur_c_cnt) {// cur_c_cnt==0, if there is no more occurence of the current character
                        if (cur_c_cnt < sigma) {
                            cur_c_cnt = cnt_c[(cur_c=alphabet[++cur_c_idx])+1];
                        }
                    }
                    size_type l=0;
                    if (i >= cnt_cc[cur_c]) { // if the current lcp entry is not already done   TODO: schleife von i bis cnt_cc[cur_c]
                        if (bwti == bwti_1 and lf < i) {  // BWT[i]==BWT[i-1]
                            l = lcp_sml[lf] ? lcp_sml[lf]-1 : 0; // l = LCP[LF[i]]-1; l < m+1
                            if (l == m) { // if LCP[LF[i]] == m+1; otherwise LCP[LF[i]] < m+1  the result is correct
                                l += (text[sai_1+m] == text[sai+m]);
                            }
                            lcp_sml[i] = l;
                            ++done_cnt;
                        } else { // BWT[i] != BWT[i-1] or LF[i] > i
                            if (lf < i)
                                l = lcp_sml[lf] ? lcp_sml[lf]-1 : 0;
                            while (text[sai_1+l] == text[sai+l] and l < m+1) {
                                ++l;
                            }
                            lcp_sml[i] = l;
                        }
                    } else { // if already done
                        l = lcp_sml[i];  // load LCP value
                    }
                    // invariant: l <= m+1
                    // begin update rmq_stack
                    size_type x = l+1;
                    size_type j = rmq_end;
                    while (x <= rmq_stack[j]) j-=2;  // pop all elements with value >= l
                    rmq_stack[++j] = i+1; // push position i
                    rmq_stack[++j] = x;	  // push value	l
                    rmq_end = j;          // update index of the value of the topmost element
                    if (lf > i) {   // if LF[i] > i, we can calculate LCP[LF[i]] in constant time with rmq
                        ++done_cnt;
                        // rmq query for lcp-values in the interval I=[prev_occ_in_bwt[BWT[i]]+1..i]
                        // rmq is linear in the stack size; can also be implemented with binary search on the stack
                        size_type x_pos = prev_occ_in_bwt[bwti]+2;
                        size_type j = rmq_end-3;
                        while (x_pos <= rmq_stack[j]) j-=2;   //  search smallest value in the interval I
                        lcp_sml[lf] = rmq_stack[j+3] - (rmq_stack[j+3]==m+2); // if lcp-value equals m+1, we subtract 1
                    }
                    if (l > m) {
                        ++nn;
                    } else
                        ++omitted_c[cur_c];

                    prev_occ_in_bwt[bwti] = i;  	 // update previous position information for character BWT[i]
                    ++cnt_cc[bwti];					 // update counter and therefore the LF information
                    sai_1 = sai;					 // update SA[i-1]
                    bwti_1 = bwti;					 // update BWT[i-1]
                }
                r_sum += r; r = bwt_buf.load_next_block();
                sa_buf.load_next_block();
            }
        }
        store_to_cache(lcp_sml, "lcp_sml", config);
    }

    // phase 2: calculate lcp_big with PHI algorithm on remaining entries of LCP
    {
        int_vector<> lcp_big(0, 0, bits::hi(n-1)+1);//nn, 0, bits::hi(n-1)+1);
        {

            mm::log("lcp-init-phi-begin");
            size_type sa_n_1 = 0;  // value for SA[n-1]
            bit_vector todo(n,0);  // bit_vector todo indicates which values are > m in lcp_sml
            {
                // initialize bit_vector todo
                int_vector_file_buffer<8> lcp_sml_buf(cache_file_name("lcp_sml", config)); // load lcp_sml
                int_vector_file_buffer<> sa_buf(cache_file_name(constants::KEY_SA, config)); // load sa
                sa_buf.load_next_block();
                for (size_type i=0, r_sum=0, r = lcp_sml_buf.load_next_block(); r_sum < n;) {
                    if (r > 0)
                        sa_n_1 = sa_buf[r-1];
                    for (; i < r_sum+r; ++i) {
                        if (lcp_sml_buf[i-r_sum] > m) {
                            todo[sa_buf[i-r_sum]] = 1;
                        }
                    }
                    r_sum += r; r = lcp_sml_buf.load_next_block();
                    sa_buf.load_next_block();
                }
            }
            rank_support_v<> todo_rank(&todo); // initialize rank for todo

            const size_type bot = sa_n_1;
            int_vector<> phi(nn, bot, bits::hi(n-1)+1); // phi

            int_vector_file_buffer<8> bwt_buf(cache_file_name(constants::KEY_BWT, config)); // load BWT
            int_vector_file_buffer<> sa_buf(cache_file_name(constants::KEY_SA, config)); // load sa
            int_vector_file_buffer<8> lcp_sml_buf(cache_file_name("lcp_sml", config)); // load lcp_sml
            sa_buf.load_next_block();
            lcp_sml_buf.load_next_block();
            uint8_t b_1 = 0;
            for (size_type i=0, r_sum=0, r=bwt_buf.load_next_block(),sai_1=0; r_sum < n;) { // initialize phi
                for (; i < r_sum+r; ++i) {
                    uint8_t b = bwt_buf[i - r_sum];  // store BWT[i]
                    size_type sai = sa_buf[i - r_sum];
                    if (lcp_sml_buf[i-r_sum] > m and b != b_1) {  // if i is a big irreducable value
                        phi[todo_rank(sai)] = sai_1;
                    } // otherwise phi is equal to bot
                    b_1 = b;
                    sai_1 = sai;
                }
                r_sum += r;
                r = bwt_buf.load_next_block();
                sa_buf.load_next_block();
                lcp_sml_buf.load_next_block();
            }
            mm::log("lcp-init-phi-end");

            mm::log("lcp-calc-plcp-begin");
            for (size_type i=0, ii=0, l=m+1,p=0; i < n and ii<nn; ++i) { // execute compact Phi algorithm
                if (todo[i]) {
                    if (i > 0 and todo[i-1])
                        l = l-1;
                    else
                        l = m+1;
                    if ((p=phi[ii]) != bot) {
                        while (text[i+l] == text[p+l]) ++l;
                    }
                    phi[ii++] = l;
                }
            }
            mm::log("lcp-calc-plcp-end");
            util::clear(text);

            mm::log("lcp-calc-lcp-begin");
            lcp_big.resize(nn);
            sa_buf.reset();
            lcp_sml_buf.reset();
            lcp_sml_buf.load_next_block();
            for (size_type i = 0, ii = 0, r = sa_buf.load_next_block(), r_sum = 0; r_sum < n and ii<nn;) {
                for (; i < r_sum+r; ++i) {
                    if (lcp_sml_buf[i - r_sum] > m) {
                        lcp_big[ii++] = phi[todo_rank(sa_buf[i - r_sum])];
                    }
                }
                r_sum += r; r = sa_buf.load_next_block();
                lcp_sml_buf.load_next_block();
            }
            mm::log("lcp-calc-lcp-end");
        }
        store_to_cache(lcp_big, "lcp_big", config);
    } // end phase 2

//		std::cout<<"# merge lcp_sml and lcp_big"<<std::endl;
    // phase 3: merge lcp_sml and lcp_big and save to disk
    {
        const size_type buffer_size = 1000000; // buffer_size has to be a multiple of 8!
        int_vector_file_buffer<> lcp_big_buf(cache_file_name("lcp_big", config)); 				// file buffer containing the big LCP values
        int_vector_file_buffer<8> lcp_sml_buf(cache_file_name("lcp_sml", config), buffer_size);// file buffer containing the small LCP values

        osfstream lcp_out_buf(cache_file_name(constants::KEY_LCP, config), std::ios::binary | std::ios::trunc);    // open out file stream
        int_vector<> lcp_buf(buffer_size, 0, lcp_big_buf.width); // buffer for the resulting LCP array
        size_type bit_size = n*lcp_big_buf.width;
        lcp_out_buf.write((char*) &(bit_size), sizeof(lcp_big_buf.int_vector_size));		   // write size
        lcp_out_buf.write((char*) &(lcp_big_buf.width), sizeof(lcp_big_buf.width));   // write int_width
        size_type wb = 0; // written bytes of data

        for (size_type i=0, r_sum=0, r=0, i2=0 ,r2_sum=0, r2=0; r_sum < n;) {
            for (; i < r_sum + r; ++i) {
                size_type l = lcp_sml_buf[i-r_sum];
                if (l > m) { // if l > m it is stored in lcp_big
                    if (i2 >= r2_sum + r2) {
                        r2_sum += r2; r2 = lcp_big_buf.load_next_block();
                    }
                    l = lcp_big_buf[i2-r2_sum];
                    ++i2;
                }
                lcp_buf[i-r_sum] = l;
            }
            if (r > 0) {
                size_type cur_wb = r*lcp_buf.width();
                if (cur_wb%8 != 0)
                    cur_wb += 7;
                cur_wb /= 8;
                lcp_out_buf.write((const char*)lcp_buf.data(), cur_wb);
                wb += cur_wb;
            }
            r_sum += r; r = lcp_sml_buf.load_next_block();
        }
        if (wb%8) {
            lcp_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
        }
    }
    sdsl::remove(cache_file_name("lcp_sml", config));
    sdsl::remove(cache_file_name("lcp_big", config));
    register_cache_file(constants::KEY_LCP, config);
    return;
}

void construct_lcp_bwt_based(cache_config& config)
{
    typedef int_vector<>::size_type size_type;
    std::string lcp_file = cache_file_name(constants::KEY_LCP, config);

    // create WT
    mm::log("lcp-bwt-create-wt-huff-begin");
    wt_huff<bit_vector, rank_support_v<>, select_support_scan<1>, select_support_scan<0>> wt_bwt;
    construct(wt_bwt, cache_file_name(constants::KEY_BWT, config));
    uint64_t n = wt_bwt.size();
    mm::log("lcp-bwt-create-wt-huff-end");

    // init
    mm::log("lcp-bwt-init-begin");
    size_type lcp_value = 0;							// current LCP value
    size_type lcp_value_offset = 0;						// Largest LCP value in LCP array, that was written on disk
    size_type phase = 0;								// Count how often the LCP array was written on disk

    size_type intervals = 0;							// number of intervals which are currently stored
    size_type intervals_new = 0;						// number of new intervals

    std::queue<size_type> q;							// Queue for storing the intervals
    std::vector<bit_vector> dict(2);							// bit_vector for storing the intervals
    size_type source = 0, target = 1;					// Defines which bit_vector is source and which is target
    bool queue_used = true;
    size_type use_queue_and_wt = n/2048;				// if intervals < use_queue_and_wt, then use queue and wavelet tree
    // else use dictionary and wavelet tree

    size_type quantity;									// quantity of characters in interval
    std::vector<unsigned char> cs(wt_bwt.sigma);				// list of characters in the interval
    std::vector<size_type> rank_c_i(wt_bwt.sigma);			// number of occurrence of character in [0 .. i-1]
    std::vector<size_type> rank_c_j(wt_bwt.sigma);			// number of occurrence of character in [0 .. j-1]

    // Calculate how many bit are for each lcp value available, to limit the memory usage to 20n bit = 2,5n byte, use at moste 8 bit
    size_type bb = (n*20-size_in_bytes(wt_bwt)*8*1.25-5*n)/n; 	// 20n - size of wavelet tree * 1.25 for rank support - 5n for bit arrays - n for index_done array
    if (n*20 < size_in_bytes(wt_bwt)*8*1.25+5*n) {
        bb = 6;
    }
    bb = std::min(bb, (size_type)8);

    size_type lcp_value_max = (1ULL<<bb)-1;		// Largest LCP value that can be stored into the LCP array
    size_type space_in_bit_for_lcp = n*bb;	// Space for the LCP array in main memory

#ifdef STUDY_INFORMATIONS
    std::cout << "# l=" << n << " b=" << (int)bb << " lcp_value_max=" << lcp_value_max << " size_in_bytes(wt_bwt<v,bs,bs>)=" << size_in_bytes(wt_bwt) << std::endl;
#endif

    // init partial_lcp
    int_vector<> partial_lcp(n, 0, bb);				// LCP array

    // init index_done
    bit_vector index_done(n, false);			// bit_vector indicates if entry LCP[i] is amend
    rank_support_v<> ds_rank_support;			// Rank support for bit_vector index_done

    // create C-array
    std::vector<size_type> C;				// C-Array: C[i] = number of occurrences of characters < i in the input
    create_C_array(C, wt_bwt);
    mm::log("lcp-bwt-init-begin-end");
    // calculate lcp
    mm::log("lcp-bwt-calc-values-begin");

    // store root interval
    q.push(0);
    q.push(n);
    intervals = 1;
    partial_lcp[0] = 0;
    index_done[0] = true;

    // calculate LCP values phase by phase
    while (intervals) {
        if (intervals < use_queue_and_wt && !queue_used) {
            mm::log("lcp-bwt-bitvec2queue-begin");
            util::clear(dict[target]);

            // copy from bitvector to queue
            size_type a2 = util::next_bit(dict[source], 0);
            size_type b2 = util::next_bit(dict[source], a2+1);
            while (b2 < dict[source].size()) {
                q.push((a2-1)>>1); q.push(b2>>1);
                // get next interval
                a2 = util::next_bit(dict[source], b2+1);
                b2 = util::next_bit(dict[source], a2+1);
            }
            util::clear(dict[source]);
            mm::log("lcp-bwt-bitvec2queue-end");
        }
        if (intervals >= use_queue_and_wt && queue_used) {
            mm::log("lcp-bwt-queue2bitvec-begin");
            dict[source].resize(2*(n+1));

            util::set_to_value(dict[source], 0);
            // copy from queue to bitvector
            while (!q.empty()) {
                dict[source][(q.front()<<1)+1 ] = 1; q.pop();
                dict[source][(q.front()<<1)   ] = 1; q.pop();
            }
            dict[target].resize(2*(n+1));

            util::set_to_value(dict[target], 0);
            mm::log("lcp-bwt-queue2bitvec-end");
        }

        if (intervals < use_queue_and_wt) {
            queue_used = true;
            intervals_new = 0;
            while (intervals) {
                // get next interval
                size_type a = q.front(); q.pop();
                size_type b = q.front(); q.pop();
                --intervals;

                wt_bwt.interval_symbols(a, b, quantity, cs, rank_c_i, rank_c_j);
                for (size_type i=0; i<quantity; ++i) {
                    unsigned char c = cs[i];
                    size_type a_new = C[c] + rank_c_i[i];
                    size_type b_new = C[c] + rank_c_j[i];

                    // Save LCP value if not seen before
                    if (!index_done[b_new] and phase == 0) {
                        partial_lcp[b_new] = lcp_value;
                        index_done[b_new] = true;
                        // Save interval
                        q.push(a_new); q.push(b_new);
                        ++intervals_new;
                    } else if (!index_done[b_new]) {
                        size_type insert_pos = b_new-ds_rank_support.rank(b_new);
                        if (!partial_lcp[insert_pos]) {
                            partial_lcp[insert_pos] = lcp_value-lcp_value_offset;
                            // Save interval
                            q.push(a_new); q.push(b_new);
                            ++intervals_new;
                        }
                    }
                }
            }
            intervals = intervals_new;
        } else {
            queue_used = false;
            intervals = 0;

            // get next interval
            size_type a2 = util::next_bit(dict[source], 0);
            size_type b2 = util::next_bit(dict[source], a2+1);

            while (b2 < dict[source].size()) {
                wt_bwt.interval_symbols(((a2-1)>>1), (b2>>1), quantity, cs, rank_c_i, rank_c_j);
                for (size_type i=0; i<quantity; ++i) {
                    unsigned char c = cs[i];
                    size_type a_new = C[c] + rank_c_i[i];
                    size_type b_new = C[c] + rank_c_j[i];
                    // Save LCP value if not seen before
                    if (!index_done[b_new] and phase == 0) {
                        partial_lcp[b_new] = lcp_value;
                        index_done[b_new] = true;
                        // Save interval
                        dict[target][(a_new<<1)+1] = 1;
                        dict[target][(b_new<<1)  ] = 1;
                        ++intervals;
                    } else if (!index_done[b_new]) {
                        size_type insert_pos = b_new-ds_rank_support.rank(b_new);
                        if (!partial_lcp[insert_pos]) {
                            partial_lcp[insert_pos] = lcp_value-lcp_value_offset;
                            // Save interval
                            dict[target][(a_new<<1)+1] = 1;
                            dict[target][(b_new<<1)  ] = 1;
                            ++intervals;
                        }
                    }
                }
                // get next interval
                a2 = util::next_bit(dict[source], b2+1);
                b2 = util::next_bit(dict[source], a2+1);
            }
            std::swap(source, target);
            util::set_to_value(dict[target], 0);
        }
        ++lcp_value;
        if (lcp_value>=lcp_value_max) {
            mm::log("lcp-bwt-write-to-file-begin");
            if (phase) {
                insert_lcp_values(partial_lcp, index_done, lcp_file, lcp_value, lcp_value_offset);
            } else {
                store_to_file(partial_lcp, lcp_file);
            }
            mm::log("lcp-bwt-write-to-file-end");
            mm::log("lcp-bwt-resize-variables-begin");
            util::init_support(ds_rank_support, &index_done); // Create rank support

            // Recalculate lcp_value_max and resize partial_lcp
            lcp_value_offset = lcp_value_max-1;
            size_type remaining_lcp_values = index_done.size()-ds_rank_support.rank(index_done.size());

            uint8_t int_width_new = std::max(space_in_bit_for_lcp / remaining_lcp_values , (size_type)bits::hi(n-1)+1);
            lcp_value_max = lcp_value_offset + (1ULL<<int_width_new);
#ifdef STUDY_INFORMATIONS
            std::cout << "# l=" << remaining_lcp_values << " b=" << (int)int_width_new << " lcp_value_max=" << lcp_value_max << std::endl;
#endif
            partial_lcp.width(int_width_new);
            partial_lcp.resize(remaining_lcp_values);
            util::set_to_value(partial_lcp, 0);
            ++phase;
            mm::log("lcp-bwt-resize-variables-end");
        }
    }
    mm::log("lcp-bwt-calc-values-end");

    // merge to file
    mm::log("lcp-bwt-merge-to-file-begin");
    if (phase) {
        insert_lcp_values(partial_lcp, index_done, lcp_file, lcp_value, lcp_value_offset);
    } else {
        store_to_file(partial_lcp, lcp_file);
    }
    register_cache_file(constants::KEY_LCP, config);
    mm::log("lcp-bwt-merge-to-file-end");
}

void construct_lcp_bwt_based2(cache_config& config)
{
    typedef int_vector<>::size_type size_type;

    uint64_t n;                         // Input length
    size_type buffer_size=1000000;      // Size of the buffer
    size_type lcp_value = 0;            // current LCP value
    std::string tmp_lcp_file = cache_file_name(constants::KEY_LCP, config)+"_tmp";
// (1) Calculate LCP-Positions-Array: For each lcp_value (in ascending order) all its occurrences (in any order) in the lcp array
    {
        mm::log("lcp-bwt2-create-wt-huff-begin");
        wt_huff<bit_vector, rank_support_v<>, select_support_scan<1>, select_support_scan<0>> wt_bwt;
        construct(wt_bwt, cache_file_name(constants::KEY_BWT, config));
        n = wt_bwt.size();
        mm::log("lcp-bwt2-create-wt-huff-begin");

        // Declare needed variables
        mm::log("lcp-bwt2-init-begin");
        size_type intervals = 0;                                        // Number of intervals which are currently stored
        size_type intervals_new = 0;                                    // Number of new intervals

        std::queue<size_type> q;                                        // Queue for storing the intervals
        std::vector<bit_vector> dict(2);                                     // bit_vector for storing the intervals
        size_type source = 0, target = 1;                               // Defines which bit_vector is source and which is target
        bool queue_used = true;                                         // Defines whether a queue (true) or the bit_vectors (false) was used to store intervals
        size_type use_queue_and_wt = n/2048;                            // if intervals < use_queue_and_wt, then use queue and wavelet tree
        // else use dictionary and wavelet tree

        size_type quantity;                                             // quantity of characters in interval
        std::vector<unsigned char> cs(wt_bwt.sigma);                         // list of characters in the interval
        std::vector<size_type> rank_c_i(wt_bwt.sigma);                       // number of occurrence of character in [0 .. i-1]
        std::vector<size_type> rank_c_j(wt_bwt.sigma);                       // number of occurrence of character in [0 .. j-1]

        // External storage of LCP-Positions-Array
        bool new_lcp_value = 0;
        uint8_t int_width = bits::hi(n)+2;

        size_type bit_size = (n+1)*int_width;                           // Size of output file in bit
        size_type wb = 0;                                               // Number of bits already written
        osfstream lcp_positions(tmp_lcp_file, std::ios::binary | std::ios::trunc);
        lcp_positions.write((char*) &(bit_size), sizeof(bit_size));     // Write length of vector
        lcp_positions.write((char*) &(int_width), sizeof(int_width));   // Write int-width of vector

        int_vector<> lcp_positions_buf(buffer_size, 0, int_width);      // Create buffer for positions of LCP entries
        size_type idx_out_buf = 0;

        bit_vector index_done(n, 0);                                    // Bitvector which is true, if corresponding LCP value was already calculated

        // Create C-array
        std::vector<size_type> C;          // C-Array: C[i] = number of occurrences of characters < i in the input
        create_C_array(C, wt_bwt);
        mm::log("lcp-bwt2-init-end");
        // Calculate LCP-Positions-Array
        mm::log("lcp-bwt2-calc-values-begin");

        // Save position of first LCP-value
        lcp_positions_buf[idx_out_buf++] = 0;
        if (new_lcp_value) {
            lcp_positions_buf[idx_out_buf-1] = lcp_positions_buf[idx_out_buf-1] + n;
            new_lcp_value = false;
        }
        if (idx_out_buf>=lcp_positions_buf.size()) {
            // Write all values from buffer to disk
            size_type cur_wb = (idx_out_buf*lcp_positions_buf.width()+7)/8;
            lcp_positions.write((const char*)lcp_positions_buf.data(), cur_wb);
            wb += cur_wb;
            idx_out_buf = 0;
        }
        index_done[0] = true;

        // Save first interval
        q.push(0); q.push(n);
        intervals = 1;

        // Calculate LCP positions
        while (intervals) {
            if (intervals < use_queue_and_wt && !queue_used) {
                mm::log("lcp-bwt2-bitvec2queue-begin");
                util::clear(dict[target]);

                // Copy from bitvector to queue
                size_type a2 = util::next_bit(dict[source], 0);
                size_type b2 = util::next_bit(dict[source], a2+1);
                while (b2 < dict[source].size()) {
                    q.push((a2-1)>>1);
                    q.push(b2>>1);
                    // Get next interval
                    a2 = util::next_bit(dict[source], b2+1);
                    b2 = util::next_bit(dict[source], a2+1);
                }
                util::clear(dict[source]);
                mm::log("lcp-bwt2-bitvec2queue-end");
            }
            if (intervals >= use_queue_and_wt && queue_used) {
                mm::log("lcp-bwt2-queue2bitvec-begin");
                dict[source].resize(2*(n+1));
                util::set_to_value(dict[source], 0);
                // Copy from queue to bitvector
                while (!q.empty()) {
                    dict[source][(q.front()<<1)+1 ] = 1; q.pop();
                    dict[source][(q.front()<<1)   ] = 1; q.pop();
                }
                dict[target].resize(2*(n+1));
                util::set_to_value(dict[target], 0);
                mm::log("lcp-bwt2-queue2bitvec-end");
            }

            if (intervals < use_queue_and_wt) {
                queue_used = true;
                intervals_new = 0;
                while (intervals) {
                    // Get next interval
                    size_type a = q.front(); q.pop();
                    size_type b = q.front(); q.pop();
                    --intervals;

                    wt_bwt.interval_symbols(a, b, quantity, cs, rank_c_i, rank_c_j);
                    for (size_type i=0; i<quantity; ++i) {
                        unsigned char c = cs[i];
                        size_type a_new = C[c] + rank_c_i[i];
                        size_type b_new = C[c] + rank_c_j[i];

                        // Save LCP value and corresponding interval if not seen before
                        if (!index_done[b_new]) {
                            // Save position of LCP-value
                            lcp_positions_buf[idx_out_buf++] = b_new;
                            if (new_lcp_value) {
                                // Mark new LCP-value
                                lcp_positions_buf[idx_out_buf-1] = lcp_positions_buf[idx_out_buf-1]+n;
                                new_lcp_value = false;
                            }
                            if (idx_out_buf>=lcp_positions_buf.size()) { // Write all values from buffer to disk
                                size_type cur_wb = (idx_out_buf*lcp_positions_buf.width()+7)/8;
                                lcp_positions.write((const char*)lcp_positions_buf.data(), cur_wb);
                                wb += cur_wb;
                                idx_out_buf = 0;
                            }
                            index_done[b_new] = true;

                            // Save interval
                            q.push(a_new);
                            q.push(b_new);
                            ++intervals_new;
                        }
                    }
                }
                intervals = intervals_new;
            } else {
                queue_used = false;
                intervals = 0;
                // get next interval
                size_type a2 = util::next_bit(dict[source], 0);
                size_type b2 = util::next_bit(dict[source], a2+1);

                while (b2 < dict[source].size()) {
                    wt_bwt.interval_symbols(((a2-1)>>1), (b2>>1), quantity, cs, rank_c_i, rank_c_j);
                    for (size_type i=0; i<quantity; ++i) {
                        unsigned char c = cs[i];
                        size_type a_new = C[c] + rank_c_i[i];
                        size_type b_new = C[c] + rank_c_j[i];
                        // Save LCP value if not seen before
                        if (!index_done[b_new]) {
                            // Save position of LCP-value
                            lcp_positions_buf[idx_out_buf++] = b_new;
                            if (new_lcp_value) {
                                // Mark new LCP-value
                                lcp_positions_buf[idx_out_buf-1] = lcp_positions_buf[idx_out_buf-1]+n;
                                new_lcp_value = false;
                            }
                            if (idx_out_buf >= lcp_positions_buf.size()) {
                                // Write all values from buffer to disk
                                size_type cur_wb = (idx_out_buf*lcp_positions_buf.width()+7)/8;
                                lcp_positions.write((const char*)lcp_positions_buf.data(), cur_wb);
                                wb += cur_wb;
                                idx_out_buf = 0;
                            }
                            index_done[b_new] = true;
                            // Save interval
                            dict[target][(a_new<<1)+1] = 1;
                            dict[target][(b_new<<1)  ] = 1;
                            ++intervals;
                        }
                    }
                    // get next interval
                    a2 = util::next_bit(dict[source], b2+1);
                    b2 = util::next_bit(dict[source], a2+1);
                }
                std::swap(source, target);
                util::set_to_value(dict[target], 0);
            }
            ++lcp_value;
            new_lcp_value = true;
        }
        mm::log("lcp-bwt2-calc-values-end");

        // Write remaining values from buffer to disk
        size_type cur_wb = (idx_out_buf*lcp_positions_buf.width()+7)/8;
        lcp_positions.write((const char*)lcp_positions_buf.data(), cur_wb);
        wb += cur_wb;
        if (wb%8) {
            lcp_positions.write("\0\0\0\0\0\0\0\0", 8-wb%8);
        }
        lcp_positions.close();
    }
// (2) Insert LCP entires into LCP array
    {
        mm::log("lcp-bwt2-reordering-begin");

        int_vector_file_buffer<> lcp_positions(tmp_lcp_file, buffer_size);

        uint8_t int_width = bits::hi(lcp_value+1)+1;             // How many bits are needed for one lcp_value?

        // Algorithm does r=ceil(int_width/8) runs over LCP-Positions-Array.
        // So in each run k>=(n/r) values of the lcp array must be calculated.
        // Because k has to be a multiple of 8, we choose number_of_values = (k+16) - ((k+16)%8)
        size_type number_of_values = ((n / ((int_width-1ULL)/8 + 1) + 16) & (~(0x7ULL)));
        int_vector<> lcp_array_out_buf(number_of_values, 0, int_width);           // Create Output Buffer

        std::string lcp_file = cache_file_name(constants::KEY_LCP, config);
        size_type bit_size = n*int_width;                               // Length of LCP-array in bit
        osfstream lcp_array(lcp_file, std::ios::binary | std::ios::trunc | std::ios::out);
        lcp_array.write((char*) &(bit_size), sizeof(bit_size));         // Write length of vector
        lcp_array.write((char*) &(int_width), sizeof(int_width));       // Write int-width of vector

        size_type wb = 0;
        for (size_type position_begin=0, position_end = number_of_values; position_begin<n and number_of_values>0; position_begin=position_end, position_end+=number_of_values) {
#ifdef STUDY_INFORMATIONS
            std::cout << "# number_of_values=" << number_of_values << " fill lcp_values with " << position_begin << " <= position <" << position_end << ", each lcp-value has " << (int)int_width << " bit, lcp_value_max=" << lcp_value << " n=" << n << std::endl;
#endif

            size_type lcp_value = 0, values=0;
            lcp_positions.reset();
            for (size_type i=0, r_sum=0, r=0; r_sum < n;) {
                for (; i < r_sum+r; ++i) {
                    size_type position = lcp_positions[i-r_sum];
                    if (position>n) {
                        position -= n;
                        ++lcp_value;
                    }
                    if (position_begin <= position and position < position_end) {
                        lcp_array_out_buf[position-position_begin] = lcp_value;
                        ++values;
                    }
                }
                r_sum += r;
                r = lcp_positions.load_next_block();
            }
            // Write next values from buffer to disk
            if (values>0) {
                size_type cur_wb = (values*lcp_array_out_buf.width()+7)/8;
                lcp_array.write((const char*)lcp_array_out_buf.data(), cur_wb);
                wb += cur_wb;
            }
        }
        // Close file
        if (wb%8) {
            lcp_array.write("\0\0\0\0\0\0\0\0", 8-wb%8);
        }
        lcp_array.close();
        register_cache_file(constants::KEY_LCP, config);
        sdsl::remove(tmp_lcp_file);
        mm::log("lcp-bwt2-reordering-end");
    } // End of phase 2
}

//void check_lcp(std::string lcpI, std::string lcpII, std::string id)
//{
//    typedef int_vector<>::size_type size_type;
//    int_vector<> lcp1,lcp2;
//    load_from_file(lcp1, (lcpI+"_"+id));
//    load_from_file(lcp2, (lcpII+"_"+id));
//    if (lcp1 != lcp2) {
//        std::cout<<"lcp results of "<<  lcpI << "and " <<lcpII<<" differ"<<std::endl;
//        for (size_type i=0, cnt=0; i<lcp1.size() and cnt<10; ++i) {
//            if (lcp1[i] != lcp2[i]) {
//                std::cout<<"i="<<i<<" "<<lcpI<<"[i]="<<lcp1[i]<<" "<<lcpII<<"[i]="<<lcp2[i]<<std::endl;
//                ++cnt;
//            }
//        }
//    }
//}
//
} // end namespace sdsl
