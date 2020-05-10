#include "xorshift.hpp"

#include <bitset>
#include <string>
#include <vector>
#include <iostream>
#include <bits/stdc++.h>

//Constructor
xorshift128Plus::xorshift128Plus(std::bitset<SEQUENCE_SIZE> seed1, std::bitset<SEQUENCE_SIZE> seed2) : m_seed1(seed1), m_seed2(seed2){};

void xorshift128Plus::seed(std::bitset<SEQUENCE_SIZE> seed1, std::bitset<SEQUENCE_SIZE> seed2)
{
    m_seed1 = seed1;
    m_seed2 = seed2;
}

std::string xorshift128Plus::xorshift(unsigned long z)
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


std::vector<std::string> xorshift128Plus::truncation(std::bitset<SEQUENCE_SIZE>& sequence, std::bitset<SEQUENCE_SIZE>& truncate,
std::bitset<BPS>& symbol, int R, int L, int MSB,std::vector<std::string>& stockage, int iteration)
{
    truncate = sequence;
    truncate >>= R;
    truncate <<= MSB - L+R;
    truncate >>= MSB - L+R;
    symbol = (truncate.to_ulong());
    std::cout<<"symboles : "<<iteration<<" || "<<symbol<<std::endl;
    stockage[iteration] = symbol.to_string();
    return stockage;
}


void xorshift128Plus::display(std::vector<std::string>& stockage,int iteration)
{
        std::cout<<" Tableau de stockage "<<iteration<<" : "<<stockage[iteration] <<std::endl;
}


// int  xorshift128Plus::channel_map(std::vector<std::string>& stockage, int iteration)
int xorshift128Plus::channel_map(std::vector<std::string>& stockage, int iteration)
{
    std::map<std::string,int> frequency;
    std::bitset<BPS> iterator(0);   

    for(int i(0); i<NB_FREQUENCIES+1;i++){   
        std::bitset<BPS> iterator(i);   
    frequency[iterator.to_string()] = 100 + i*100;}
    	std::cout<<"frequence [964] = "<<frequency["1110001"]<<std::endl;
    return frequency[stockage[iteration]];
}
