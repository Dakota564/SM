#include "xoroshiro256plusplus.hpp"

#include <bitset>
#include <string>
#include <vector>
#include <iostream>
#include <bits/stdc++.h>

static uint64_t s[4]={1,1,1,1};
// static uint64_t s[4]={1,0,0,0};

//Constructor
xoroshiro256plusplus::xoroshiro256plusplus(uint32_t test, std::bitset<SEQUENCE_SIZE> seed1, std::bitset<SEQUENCE_SIZE> seed2) : test(0xfeedface00000000ULL | test),m_seed1(seed1), m_seed2(seed2){};

// void xoroshiro256plusplus::seed(std::bitset<SEQUENCE_SIZE> seed1, std::bitset<SEQUENCE_SIZE> seed2)
// {
//     m_seed1 = seed1;
//     m_seed2 = seed2;
// }s

uint64_t *xoroshiro256plusplus::seed(uint64_t *array)
{
        int value =31;  // assuming a 32 bit int
        int p;
        int k;
		int cnt0;
		int cnt1;
        //  static uint64_t array[64] = {0};

          for (p = 0,k=31; p <32,k>=0 ;++p,--k) {
			  array[k] = (value >> p) & 1;
			  std::cout<<array[p];
			  if (array[p]==0){
				 cnt0++;			  
			  }
			  else 
			  	cnt1++;
				}
		std::cout<<"number of 0 arr: "<<cnt0<<std::endl;
		std::cout<<"number of 1 arr: "<<cnt1<<std::endl;
		 return array;
}



uint32_t xoroshiro256plusplus::lehmer64(void) {
  test *= 2685821657736338717ULL;
  std::cout<<"Lehmer function : "<<(test >> 32)<<std::endl;
  return uint32_t(test >> 32);
}




// static uint64_t s[4]={1,0,1,0};

uint64_t xoroshiro256plusplus::next256plusplus(uint64_t *s) {
	// std::cout<<"Voici s[0] "<<s[0]<<s[1]<<s[2]<<s[3]<<std::endl;
	const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
	const uint64_t t = s[1] << 17;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = rotl(s[3], 45);
	// std::cout<<"Voici le result :"<<result<<std::endl;
	return result;
}

uint64_t xoroshiro256plusplus::next128starstar(uint64_t *s) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = rotl(s0 * 5, 7) * 9;

	s1 ^= s0;
	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotl(s1, 37); // c

	return result;
}

uint64_t xoroshiro256plusplus::next128plusplus(uint64_t *s) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = rotl(s0 + s1, 17) + s0;

	s1 ^= s0;
	s[0] = rotl(s0, 49) ^ s1 ^ (s1 << 21); // a, b
	s[1] = rotl(s1, 28); // c

	return result;
}

uint64_t xoroshiro256plusplus::next256starstar(uint64_t *s) {
	const uint64_t result = rotl(s[1] * 5, 7) * 9;

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void xoroshiro256plusplus::jump(void) {
	static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= *s;
				s1 ^= *(s+1);
				s2 ^= *(s+2);
				s3 ^= *(s+3);
			}
			// next(s);	
			// next256starstar(s);	
			next128starstar(s);	
			
		}
		
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}


void xoroshiro256plusplus::long_jump(void) {
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			// next(s);	
		}
		
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}






std::string xoroshiro256plusplus::xorshift(unsigned long z)
{
    while (z--){
    std::bitset<SEQUENCE_SIZE> t = m_seed1;
    m_seed1 = m_seed2;
    t ^= t << 23;
    t ^= t >> 17;
    t ^= m_seed2;
    t ^= m_seed2 ^ (m_seed2 >> 26);
    m_seed1 = t;}
    return (m_seed1.to_string() + m_seed2.to_string());}

























std::vector<std::string> xoroshiro256plusplus::truncation(std::bitset<SEQUENCE_SIZE>& sequence, std::bitset<SEQUENCE_SIZE>& truncate,
std::bitset<BPS>& symbol, int R, int L, int MSB,std::vector<std::string>& stockage, int iteration)
{
    truncate = sequence;
    truncate >>= R;
    truncate <<= MSB - L+R;
    truncate >>= MSB - L+R;
    symbol = (truncate.to_ulong());
    stockage[iteration] = symbol.to_string();
	// std::cout<<"STOCKAGE : "<<iteration<<" == "<<stockage[iteration]<<std::endl;
    return stockage;
}


void xoroshiro256plusplus::display(std::vector<std::string>& stockage)
{
 for (int i(0); i<130; i++){
        std::cout<<" Tableau de stockage "<<i<<" : "<<stockage[i] <<std::endl;
        }
}

int xoroshiro256plusplus::channel_map(std::vector<std::string>& stockage, int iteration)
{
    std::map<std::string,int> frequency;
    std::bitset<BPS> iterator(0);   

    for(int i(0); i<NB_FREQUENCIES+1;i++){   
        std::bitset<BPS> iterator(i);   
    	frequency[iterator.to_string()] = 100 + i*100;}
    return frequency[stockage[iteration]];
}
