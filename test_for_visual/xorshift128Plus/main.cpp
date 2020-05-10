#include "xorshift.cpp"
#include <iostream>
#include <bitset>


int main(int argc, char const *argv[])
{
    xorshift128Plus x128;

    x128.seed(seed1,seed2);                                        //Initialisation of the seeds
    std::bitset<SEQUENCE_SIZE> sequence(x128.xorshift(NB_OPERATIONS));     //Determination of the Pseudo Random sequence    
    while(iteration<NB_HOPS){                              //Operations for NB_FREQUENCIES
     //    for(loop2=1; loop2<REARRANGEMENT;loop2++){ //Xorshit of a frequency's sequence
          std::bitset<SEQUENCE_SIZE> sequence(x128.xorshift(105)); //Generate the xorshift sequence 
          x128.truncation(sequence, truncate, symbol,R,L,MSB,stockage,iteration); //Generate the binary output (the random number)
     //    }


        //Mapping binary outputs to frequencies
         stockage_frequency[iteration] = x128.channel_map(stockage,iteration);
     // x128.channel_map(stockage,iteration);
        std::cout<<"Frequency : "<<iteration<<" ||  "<<stockage[iteration]<<" || "<<stockage_frequency[iteration]<<std::endl;
     // std::cout<<"Stockage iteration = "<<stockage_frequency[iteration]<<std::endl;
        iteration++;
    }
     // ------------- TEST PART ----------------

    // int ctr[NB_FREQUENCIES]={0};
    int ctr(0);

    // for (int i=0; i<NB_FREQUENCIES; i++){
     //   std::cout<<"nb of frequency : "<<stockage_frequency[i]<< std::endl;

    //  if (stockage_frequency[i] == 100){
            //  std::cout<<"BOUAUAAUAU : "<<ctr << std::endl;
            // bool exists = std::find(std::begin(stockage_frequency), std::end(stockage_frequency), 600) != std::end(stockage_frequency);
            // bool exists = std::find(std::begin(stockage_frequency), std::end(stockage_frequency), 600) != std::end(stockage_frequency);
            // while (exists)
            //           ctr++;
//     r(*(ptr_ctr+i))++;}

    // }
    int variance[NB_FREQUENCIES]={0};
     for(int i(0);i<NB_FREQUENCIES+1;i++){
        std::cout<<i+1<<" : "<<std::count(std::begin(stockage_frequency), std::end(stockage_frequency), 100 + i*100) << std::endl;
        variance[i]= abs((NB_HOPS/NB_FREQUENCIES)-std::count(std::begin(stockage_frequency), std::end(stockage_frequency), 100 + i*100));
        
        }
     int sum(0);
     sum = std::accumulate(std::begin(variance), std::end(variance),sum);
     std::cout<<"The max variation is :"<<*std::max_element(std::begin(variance), std::end(variance))<<std::endl;
     std::cout<<"The total variation is :"<<sum<<std::endl;
     std::cout<<"FREQUENCE 964 :"<<stockage[964]<<std::endl;


    return 0;
}