#ifndef XORSHIFT256PLUSPLUS_HPP
#define XORSHIFT256PLUSPLUS_HPP

#include <bitset>
#include <string>
#include <vector>
#include <iostream>
#include <cstdint>
#include <cstdlib>

// typedef unsigned __int128 uint128_t;
uint64_t high; 
uint64_t low;
unsigned long long test((((uint64_t) high) << 32) | ((uint64_t) low));


//================ Variable declaration ===================
const int SEQUENCE_SIZE = 64;                        //Number of bits in a sequence
const int BPS           = 7; //Bits per symbol
const int seed1         = 100;                       //First seed to feed the xorshift
const int seed2         = 20;                        //Second seed to feed the xorshift
const int R             = 50;      //Right shift ||   48
const int L             = R+BPS;                     //Left shift
 
unsigned long NB_OPERATIONS = 32;  //100  Define the number of xorshift to generate a sequence  ||   256
const int REARRANGEMENT     = 2;   //Optimal number of xorshift of the sequence for almost equiprobable frequency occurence ||35
const int NB_FREQUENCIES    = 79;   //Number of frequencies generate by the PN generator
const int NB_HOPS           = 500; //Number of hops generate by the PN generator

// static uint64_t s[2]={0,0};


std::bitset<SEQUENCE_SIZE> truncate;
std::bitset<BPS> symbol;
std::vector<std::string> reduce_sequence;
std::vector<std::string> stockage(10000);
std::vector<int> stockage_frequency(10000);
int MSB(SEQUENCE_SIZE);
int iteration(0);
//=========================================================
// int x(10),k(20);
class xoroshiro256plusplus
{
  private:
    std::bitset<SEQUENCE_SIZE> m_seed1;
    std::bitset<SEQUENCE_SIZE> m_seed2;
    uint64_t  test;

  public:
    xoroshiro256plusplus(uint32_t test = 0xcafef00d,std::bitset<SEQUENCE_SIZE> seed1=0, std::bitset<SEQUENCE_SIZE> seed2=0);
    uint32_t lehmer64(void);
    uint64_t *seed(uint64_t *array);
    uint64_t rotl(const uint64_t x, int k) {return (x << k) | (x >> (64 - k));}
    // uint64_t rotl(const uint64_t x, int k);
    uint64_t next128starstar(uint64_t *s);
    uint64_t next128plusplus(uint64_t *s);
    uint64_t next256plusplus(uint64_t *s);
    uint64_t next256starstar(uint64_t *s);

    void jump(void);
    void long_jump(void);
    std::string xorshift(unsigned long z);
    std::vector<std::string> truncation(std::bitset<SEQUENCE_SIZE>& sequence, std::bitset<SEQUENCE_SIZE>& truncate, std::bitset<BPS>& symbol, int R, int L, int MSB,std::vector<std::string>& stockage, int iteration);
    void display(std::vector<std::string>& stockage);
    int channel_map(std::vector<std::string>& stockage, int iteration);
};
#endif