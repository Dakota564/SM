#include "xoroshiro256plusplus.cpp"
#include <iostream>
#include <bitset>
#include <string>
#include <typeinfo>
#include <cstdint>
#include <cstdlib>


int main(int argc, char const *argv[])
{
    xoroshiro256plusplus x128;


//     x128.seed(seed1,seed2);                                        //Initialisation of the seeds
    //std::bitset<SEQUENCE_SIZE> sequence(x128.xorshift(NB_OPERATIONS));     //Determination of the Pseudo Random sequence    
    while(iteration<= NB_HOPS){                              //Operations for NB_FREQUENCIES
      //   for(int loop2=1; loop2<REARRANGEMENT;loop2++){ //Xorshit of a frequency's sequence

         //  std::bitset<SEQUENCE_SIZE> sequence(x128.next(s)); //Generate the xorshift sequence 
         //  std::cout<<x128.lehmer64()<<std::endl; //Generate the xorshift sequence 

          

        //   std::bitset<SEQUENCE_SIZE> sequence(x128.lehmer64()); //Generate the xorshift sequence 

          // std::bitset<SEQUENCE_SIZE> sequence(x128.next128plusplus(s)); //Generate the xorshift sequence 
          // std::bitset<SEQUENCE_SIZE> sequence(x128.next128starstar(s)); //Generate the xorshift sequence 
            //  std::bitset<SEQUENCE_SIZE> sequence(x128.next256plusplus(s)); //Generate the xorshift sequence 
          std::bitset<SEQUENCE_SIZE> sequence(x128.next256starstar(s)); //Generate the xorshift sequence 


         //  std::cout<<"Voici la sequence : "<<sequence<<std::endl;
          // x128.jump();


          // std::bitset<SEQUENCE_SIZE> sequence(*s); //Generate the xorshift sequence 
          // std::cout<<"Voici le jump : "<<sequence<<std::endl;
     // //     uint64_t *s =  (x128.seed());   
                                  
          
          x128.truncation(sequence, truncate, symbol,R,L,MSB,stockage,iteration); //Generate the binary output (the random number)

     // x128.jump();
     //    //Mapping binary outputs to frequencies

        stockage_frequency[iteration] = x128.channel_map(stockage,iteration);
             std::cout<<"Frequency : "<<iteration<<" ||  "<<stockage[iteration]<<" || "<<stockage_frequency[iteration]<<std::endl;

     //    iteration++;



          
      //   }
               //     std::cout<<"Sequence : "<<stockage[iteration]<<std::endl;

 iteration++;
    }
    
    

     // ------------- TEST PART ----------------
int variance[NB_FREQUENCIES]={0};
     for(int i(0);i<NB_FREQUENCIES;i++){
        std::cout<<i+1<<" : "<<std::count(std::begin(stockage_frequency), std::end(stockage_frequency), 100 + i*100) << std::endl;
        variance[i]= abs((NB_HOPS/NB_FREQUENCIES)-std::count(std::begin(stockage_frequency), std::end(stockage_frequency), 100 + i*100));
        
        }
     int sum(0);
     sum = std::accumulate(std::begin(variance), std::end(variance),sum);
     std::cout<<"The max variation is :"<<*std::max_element(std::begin(variance), std::end(variance))<<std::endl;
     std::cout<<"The total variation is :"<<sum<<std::endl;



     
   
   //  std::cout<<sizeof(test);
    return 0;

}