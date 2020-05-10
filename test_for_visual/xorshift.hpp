#ifndef XORSHIFT_HPP
#define XORSHIFT_HPP

#include <bitset>
#include <string>
#include <vector>
#include <iostream>

//================ Variable declaration ===================
const int BPS           = 2;        //Bits per symbol
const int seed1         = 100;      //First seed to feed the xorshift
const int seed2         = 20;       //Second seed to feed the xorshift
const int R             = 48;       //Right shift
const int L             = R+BPS;    //Left shift
const int SEQUENCE_SIZE = 64;       //Number of bits in a sequence
 
unsigned long NB_OPERATIONS = 1;  //100  Define the number of xorshift to generate a sequence
const int REARRANGEMENT     = 35;   //Optimal number of xorshift of the sequence for almost equiprobable frequency occurence 
const int NB_FREQUENCIES    = 100; //Number of frequencies generate by the PN generator

std::bitset<SEQUENCE_SIZE> truncate;
std::bitset<BPS> symbol;
std::vector<std::string> reduce_sequence;
std::vector<std::string> stockage(10000);
std::vector<int> stockage_frequency(10000);
int MSB(SEQUENCE_SIZE);
int iteration(0);
//=========================================================

class xorshift128Plus
{
  private:
    std::bitset<SEQUENCE_SIZE> m_seed1;
    std::bitset<SEQUENCE_SIZE> m_seed2;
  public:
    xorshift128Plus(std::bitset<SEQUENCE_SIZE> seed1=0, std::bitset<SEQUENCE_SIZE> seed2=0);
    void seed(std::bitset<SEQUENCE_SIZE> seed1, std::bitset<SEQUENCE_SIZE> seed2);
    std::string xorshift(unsigned long z);
    std::vector<std::string> truncation(std::bitset<SEQUENCE_SIZE>& sequence, std::bitset<SEQUENCE_SIZE>& truncate, std::bitset<BPS>& symbol, int R, int L, int MSB,std::vector<std::string>& stockage, int iteration);
    void display(std::vector<std::string>& stockage);
    int channel_map(std::vector<std::string>& stockage, int iteration);
};
#endif