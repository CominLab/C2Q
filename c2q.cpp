
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include<iostream>
#include<fstream>
#include<cstring>
//#include<cmath>
#include<vector>

# include <stdio.h>
# include <stdlib.h>

#include<algorithm>
#include<ext/hash_map>
using __gnu_cxx::hash_map;
using namespace std;

enum FastLineType {LineHeader, LineBases, LineSuppHeader, LineQualities};


// In order to give more precise calculation and avoid zero in dominator, 
// we introduce the SCIENTIFIC_NUMBER calculation. 
// We will re-define addtion, multiplication and power in the following sub rountine functions.
// It is contributed by Prof. Minghua Deng, Peking University.

struct SCIENTIFIC_NUMBER
{
	int factor;
	double value;
};

int NS = 0;   
int k = 0;  
int d = 0;   // d = 0: take the whole data
int M = 0;
char type = 'A';


// Save the data file directions and names
char* fileDirect[200];
char* fileName[200];


unsigned long power;//parameter, power = 4^(k-1)

//***change para  : the Matrix for save c2, c2star, c2shepp matrix.   
//*******The maximum number of species(datasets) is 200. Users can change the parameters here.
double c2Form[200][200], c2starForm[200][200], c2sheppForm[200][200]; 
double c2QForm[200][200], c2QstarForm[200][200];
double d2Form[200][200], d2starForm[200][200], d2sheppForm[200][200];
double d2QForm[200][200], d2QstarForm[200][200];
    
// these hash tables are used as temporary tables while computing Xw for a single dataset
hash_map<unsigned long,unsigned long> HashTable;
hash_map<unsigned long,double> HashTableQ;

// these hash tables contains the kmers count for all datasets
hash_map<unsigned long,unsigned long> HashTableS[200];
hash_map<unsigned long,double> HashTableSQ[200];

//receptacle of Pw(probability of a kmer word) 
hash_map<unsigned long,SCIENTIFIC_NUMBER > HashPw;

//***change para : receptacle of Pw(probability of a kmer word) for all datasets
hash_map<unsigned long,SCIENTIFIC_NUMBER > HashPwS[200];

int seqlength=99999;

// seq: save read/genome(A C G T) in each line from file
char seq[99999];

// inverse of seq
char seq_inverse[99999];

// sequence of quality values
char qual[99999];

// this is the quality encoding offset change it for different encoding schemes
char offset = '!';


// Prior[r][m][n][l] records the probability of each kmer word(Pw) in rth dataset. 
// Each kmer word can be decomposed into m As, n Cs, l Gs and k-m-n-l Ts.
// 30 is not related to 30 species. //***change para
// *******The maximum number of k is 29. Users can change this parameter.
SCIENTIFIC_NUMBER Prior[200][30][30][30];  


// records the frequency of A, C, G, T in a dataset
double anum=0, cnum=0, gnum=0, tnum=0, nnum=0;



// Scientific Number calculation funtion:

// TransToReal : trans a SCIENTIFIC NUMBER to a read number
double TransToReal(SCIENTIFIC_NUMBER dSci)
{
    double dReal=0;
    dReal = dSci.value * pow(10,dSci.factor);
    return dReal;

}


// TransToScientific : trans a read number to a SCIENTIFIC NUMBER
SCIENTIFIC_NUMBER TransToScientific(double dReal)
{
	SCIENTIFIC_NUMBER sciTemp;
	int count;

	if( dReal==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
	}
	else if(dReal>10.0 || dReal<-10.0)
	{
		count=0;
		while(dReal>10.0 || dReal<-10.0)
		{
			dReal /=10.0;
			count++;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else if( dReal<1.0 && dReal>-1.0)
	{
		count=0;
		while( dReal<1.0 && dReal>-1.0 )
		{
			dReal *=10.0;
			count--;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else
	{
		sciTemp.value=dReal;
		sciTemp.factor=0;
	}

	return sciTemp;
}


// SciMultiple : Multiplication of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
//    cout << "SciMultiple " << endl;
    
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int count;

	if( left.value==0.0 || right.value==0.0 )
	{
//        cout << "Both 0 " << endl;
        
		sciTemp.value=0.0;
		sciTemp.factor=0;

		return sciTemp;
	}

	// now both left and right element are nonzero
	dTemp=left.value * right.value;
    
//    cout << "left.value " << left.value << endl;
//    cout << "right.value " << right.value << endl;
//    cout << "dTemp " << dTemp << endl;
    
    
	if( dTemp>10.0 || dTemp<-10.0 )
	{
        
//        cout << "10 < dTemp or dTemp < -10 " << endl;
        
		count=0;
		while(dTemp>10.0 || dTemp<-10.0 )
		{
			dTemp /=10.0;
			count++;
		}
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
	}
	else if( dTemp<1.0 && dTemp>-1.0)
	{
//        cout << "dTemp < 1 or dTemp > -1 " << dTemp << endl;
        
		count=0;
		while( dTemp<1.0 && dTemp>-1.0 )
		{
			dTemp *=10.0;
			count--;
		}
        
//        cout << "count " << count << endl;
        
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
        
//        cout << "sciTemp " << sciTemp.value << " " << sciTemp.factor << endl;
	}
	else
	{
        
//        cout << "dTemp normal " << dTemp << endl;
        
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor;
	}

	return sciTemp;
}


// SciMultiple : Multiplication between a SCIENTIFIC NUMBERS and a read number
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;

	sciTemp=TransToScientific(right);
	sciTemp=SciMultiple(left,sciTemp);

	return sciTemp;
}





// SciAddition : addition of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int i,count;

	if( left.value==0.0 || right.value==0.0 )
	{
		if( left.value==0.0 )
			return right;
		else
			return left;
	}

	// now the two element are both non zero
	if( left.factor>=right.factor)
	{
		// left element is larger than right element
		dTemp=right.value;
		for(i=0;i<(left.factor-right.factor);i++)
			dTemp /=10.0;
		dTemp +=left.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}

		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while(dTemp>10.0 || dTemp<-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor;
		}
		return sciTemp;
	}
	else
	{
		// right element  is larger than left element
		dTemp=left.value;
		for(i=0;i<(right.factor-left.factor);i++)
			dTemp /=10.0;
		dTemp +=right.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}

		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while( dTemp>10.0 || dTemp <-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor;
		}
		return sciTemp;
	}
}


// SciAddition : addition between a SCIENTIFIC NUMBERS and a real number
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;

	sciTemp=TransToScientific(right);
	sciTemp=SciAddition(left,sciTemp);

	return sciTemp;
}

// SciPow : give the power of a scientific number
SCIENTIFIC_NUMBER SciPow(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
	double dTemp;
	int iTemp;
/*
	if(left.value==0.0 )
	{
		printf("the base of the power is nagative\n");
		exit(1);
	}
*/
	if(left.value==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
		return sciTemp;
	}

	dTemp=(log10(fabs(left.value))+left.factor)*right;

	if( dTemp>0.0 )
		iTemp=int(ceil(dTemp));    //ceil(a)是求不小于a的最小整数。floor(a)表示求不大于a的最大整数
	else
		iTemp=int(floor(dTemp));
	sciTemp.factor=iTemp;
	sciTemp.value=pow(10.0,dTemp-iTemp);

	return sciTemp;
}

// given a seq computes the score of subword seq_Q[index, ..., index+k-1]
double computeQualityScore( char* seq_Q, int index ) {
  double quality = 1;
  for( int i = index; i < index + k; ++i) {
    double quality_value = (double) seq_Q[i] - (double) offset;
    // !!! Here a lookup table can be used to speed up the whole computation
    double probQ = 1 - pow(10, -quality_value/10 ); 
    quality = quality*probQ;
    //    cout << "* " << probQ << " " << quality_value << " " << seq_Q[i] <<  i << endl;

  }
  return quality;
}

double computeAWP(int dataSet, unsigned long index) {
  double den = (double)HashTableS[dataSet][index];
  double num = HashTableSQ[dataSet][index] ;
  return ( (den > 0) ? num / den : 0 );
}


// computes the reverse complement of an unsigned long representing a kmer
// note that when the bit representation of bases is {A=00, C=01, G=10, T=11}
// reverse complement of x is
//    Reverse(Not(x))
// where both Reverse and Not must be interpreted as bitwise operations
unsigned long reverseComplement(unsigned long index, int k) {
  unsigned long rc = 0;
  index = ( ~(index) & ( (1 << 2*k) - 1) );
  for (int i = 0; i < k; i++) {
    rc <<= 2; // equivalent to (rc = rc * 4)
    rc += (index & 0x3); // equivalent to (rc = index % 4)
    //rc <<= 2; // equivalent to (rc = rc * 4)
    index >>= 2; // equivalent to (index = index / 4)    
  }
  return rc;
}

// The algorithm in SeqKmerCount: to count the frequency of each kmer in a seqence of read/genome
//	k=6
//	ACGTCACGTACGT...
//	ACGTCA index1
//	 CGTCAC index2 = floor(index1/4^5)*4 + 1(C)
//	  GTCACG index3 = floor(index2/4^5)*4 + 2(G)
//	   TCACGT index4 = floor(index3/4^5)*4 + 3(T)
//	    CACGTA index5 = floor(index4/4^5)*4 + 0(A)
//	     ACGTAC index6 = ...
unsigned long SeqKmerCount(char* seq, int k, char* qual = NULL) {
  // int count: The length of char seq 
  int count = 0;     

  int i=0, j=0;
  unsigned long index = 0;
  unsigned long total = 0; //total number of the kmers counted in seq
  
  while(seq[i]) {
    //kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
    //current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
    if(seq[i]=='A'|| seq[i] == 'a') {j++; anum++;}
    else if(seq[i]=='C'|| seq[i] == 'c') { j++; index++; cnum++;}
    else if(seq[i]=='G'|| seq[i] == 'g') { j++; index+=2; gnum++;}
    else if(seq[i]=='T'|| seq[i] == 't') { j++; index+=3; tnum++;}
    else { j=0; index=0; nnum++;}//If seq[i] is ambiguous, reset j and index
    
    if( j == k ) {
      HashTable[index]++;
      if (qual != NULL) {
	double P = computeQualityScore(qual,j-k+1);
	HashTableQ[index] += P;
	HashTableQ[reverseComplement(index, k)] += P;
      }
      total++;        
      index %= power;// current index = floor(previous index/4^(k-1))
      j--;//the lengh of seq[i+1,i+2,...i+k-1]
    }
    index*=4;//current index = floor(previous index/4^(k-1))*4
    i++;
    count++;
  }
  // Inversed seq from right to left, but not inverse A-T, G-C
  for (int m = 1; m != count + 1; m++) {
    seq_inverse[m - 1] = *(seq + count - m);
  }
  seq_inverse[count] = '\0';
  i=0,j=0;
  index = 0;

  while(seq_inverse[i])	{
    // Kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
    // Current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
    if(seq_inverse[i]=='T'|| seq_inverse[i] == 't') {j++;}
    else if(seq_inverse[i]=='G'|| seq_inverse[i] == 'g') { j++; index++;}
    else if(seq_inverse[i]=='C'|| seq_inverse[i] == 'c') { j++; index+=2;}
    else if(seq_inverse[i]=='A'|| seq_inverse[i] == 'a') { j++; index+=3;}
    else { j=0; index=0; nnum++;}// If seq[i] is ambiguous, reset j and index
    
    if( j == k ) {
      HashTable[index]++;
      index %= power;// Current index = floor(previous index/4^(k-1))
      j--;// The lengh of seq[i+1,i+2,...i+k-1]
    }
    index*=4;// Current index = floor(previous index/4^(k-1))*4
    i++;
  }
  return total;      
}


// This function scans the input reads file and invokes the proper
// function (i.e. either with or without qualities) to compute the
// number of kmers with or without quality weight respectively
unsigned long DataKmerCount(int DataNum, char* argv_DataNum, int k) {

    // Initial HashTable
    HashTable.clear();   
    HashPw.clear();
    HashTableQ.clear();

    // Initialize total_DataNum. total_DataNum records the total number of kmer in rth dataset
    unsigned long total_DataNum = 0;

    // Initialize anum, cnum, gnum, tnum, recording the frequency of A, T, C, G of rth dataset
    anum=0; cnum=0; gnum=0; tnum=0;

    // Initialize seq with zeros
    memset(seq, 0, seqlength * sizeof(char));
    memset(qual, 0, seqlength * sizeof(char));
        
    // Parameters: k and power = 4^(k-1)
    //    power = 1; for( int i = 0; i < k-1; i++) power *= 4;   
    power = 1 << (2 * (k - 1) );

    // Open Data from file argv[r*2-1]
    ifstream fin(argv_DataNum);   

    // Scanning was completely wrong on the original one, the following code
    // solves the problem and take into account quality values when required
    // This solution now works only if input file is well formed according to
    // the type passed (ie A for fasta and Q for fastq)
    FastLineType lineType = LineHeader;
    char seqTmp[99999];
    while(fin.getline(seqTmp,seqlength)) 
    {
      // all this if and  else if could be substituted with a 'switch' statement
      if (lineType == LineHeader) {
	// on the header we skip over the next line that is supposed to contain
	// the bases.
	lineType = LineBases;
      }
      else if(lineType == LineBases) {	
	strcpy(seq, seqTmp);
	if ( (type == 'A') || (type == 'a') ) {
	  total_DataNum += SeqKmerCount(seq,k, NULL);
	  lineType = LineHeader;
	} else {
	  lineType = LineSuppHeader;
	}
      }
      else if(lineType == LineSuppHeader) {
	lineType = LineQualities;
      }
      else if(lineType == LineQualities) {
	total_DataNum += SeqKmerCount(seq, k, seqTmp);
	lineType = LineHeader;
      }	          
    }
    fin.close();


     // Sort the kmer
     vector<unsigned long> temp;
     for( hash_map<unsigned long,unsigned long>::iterator i = HashTable.begin(); i!= HashTable.end(); i++) temp.push_back(i->first);
     sort(temp.begin(), temp.end());


     unsigned long key;
     for(vector<unsigned long>::iterator j = temp.begin(); j!=temp.end(); j++)
     {
         key = *j;
	 HashTableS[DataNum][key]=HashTable[key];
	 HashTableSQ[DataNum][key] = HashTableQ[key];
//	 fout2.write((char*)&key,sizeof(unsigned long));
//	 fout2.write((char*)&HashTableS[DataNum][key],sizeof(unsigned long));
//	 fout2 << "Key" << key;
//	 fout2 << key << "," << HashTable[key] << endl;
//	 fout2 << "Count" << HashTableS[DataNum][key] << endl;
     }
//	fout2.close();


     return total_DataNum;

}



// Compute Pw for each kmer word
void ComputePw(int DataNum, int k, char* argg_DataNumK, char* argg_DataNum, int d)
{
     
     //Compute probability(frequency) of A, C, G, T: pa, pc, pg, pt

     // Initialize pa, pc, pg, pt 
     SCIENTIFIC_NUMBER pa, pc, pg, pt;
     pa.value = 0; pc.value = 0; pg.value = 0; pt.value = 0;
     pa.factor = 0; pc.factor = 0; pg.factor = 0; pt.factor = 0;

     // Total count of A, C, G, T
     double totalACGT = 0;
     totalACGT = anum + cnum + gnum + tnum;

     // Frequency of A, C, G, T
     pa=TransToScientific(anum/totalACGT);
     pc=TransToScientific(cnum/totalACGT);
     pg=TransToScientific(gnum/totalACGT);
     pt=TransToScientific(tnum/totalACGT);
//   cout<<pa<<pc<<pg<<pt<<endl;

 
     // Output probability of A, C, G, T
     char dstr[5]; sprintf(dstr, "%d", d); 
     char outputfile_pa[1000];
     strcpy(outputfile_pa,argg_DataNum);
     strcat(outputfile_pa,"_k");
     strcat(outputfile_pa,argg_DataNumK);
     strcat(outputfile_pa,"_d"); 
     strcat(outputfile_pa,dstr);
     strcat(outputfile_pa,"_pa"); // strcpy(outputfile,argv[1]) strcat(outputfile,argv[2])
     ofstream fout_pa(outputfile_pa);

     // Convert SCITIFIC_NUMBER to real number
     long double pa_real = TransToReal(pa);
     long double pc_real = TransToReal(pc);
     long double pg_real = TransToReal(pg);
     long double pt_real = TransToReal(pt);
     fout_pa << "pa" << pa_real << endl;
     fout_pa << "pc" << pc_real << endl;
     fout_pa << "pg" << pg_real << endl;
     fout_pa << "pt" << pt_real << endl << endl;
     fout_pa.close();


     // Computer Pw(probability) for all kmer word

     // First compute A[] C[] G[] T[]: A[pa^0, pa^1, pa^2, ... , pa^k] 
     SCIENTIFIC_NUMBER A[k+1], C[k+1], G[k+1], T[k+1];
     for (int p = 0; p < k + 1; p++)
     {
//         cout << "p" << p << " k" << k << endl;
         
         A[p]=SciPow(pa,p);
//         cout<< "A[p]" << A[p].value << A[p].factor << endl;
         C[p]=SciPow(pc,p);
//         cout<< "C[p]" << C[p].value << C[p].factor << endl;
         G[p]=SciPow(pg,p);
//         cout<< "G[p]" << G[p].value << G[p].factor << endl;
         T[p]=SciPow(pt,p);
//         cout<< "T[p]" << T[p].value << T[p].factor << endl;
     }
//     cout<<endl;


     // Second calculate probability for each kmer word: 
     // each kmer word can be decomposed into n_1 As, n_2 Cs, n_3 Gs, and n_4 Ts, such that n_1 + n_2 + n_3 + n_4 = k
     int m, n, l;
     for( m = 0; m <k+1; m++)
     {
      //   Prior[m][n][l] = A[m];
         for( n = 0; n < k+1 - m; n++)
         {
          //   Prior[m][n][l]*=C[n];
             for( l = 0; l < k+1 -m - n; l++)
             {
                    Prior[DataNum][m][n][l] = SciMultiple(SciMultiple(SciMultiple(A[m], C[n]), G[l]), T[k-m-n-l]);
		    //                    cout << "m" << m << "n" << n << "l" << l << " Prior " <<  Prior[DataNum][m][n][l].value << " " << Prior[DataNum][m][n][l].factor <<endl;
//                    cout << "Prior-020 " << Prior[DataNum][0][2][0].value << " " << Prior[DataNum][0][2][0].factor <<endl;
             }
         }
     }
//   cout<<endl;


}


// Print kmer count hashtable and its corresponding Pw
void PrintKmerCountPw(int k, int DataNum, char* argg_DataNumK, char* argg_DataNum, int d)
{
     
     // Output kmer count-pw files
     char dstr[5]; sprintf(dstr, "%d", d); 
     char outputfile_CountPw[1000]; 
     strcpy(outputfile_CountPw,argg_DataNum); 
     strcat(outputfile_CountPw,"_k");
     strcat(outputfile_CountPw,argg_DataNumK);
     strcat(outputfile_CountPw,"_d"); // strcpy(outputfile,argv[1]) strcat(outputfile,argv[2])
     strcat(outputfile_CountPw,dstr);
     strcat(outputfile_CountPw,"_wordcount_pw"); // strcpy(outputfile,argv[1]) strcat(outputfile,argv[2])
     ofstream fout_CountPw(outputfile_CountPw);


     for(unsigned long kk = 0; kk < pow(4,k); kk++)
     {

         // For each kmer, convert base-10 numeral system back to base-4 numeral system
         int a; int b[30]; memset(b, 0, 30);
         int t=k-1;
         a = kk; 
         while(a>=3) {b[t]=a%4;a/=4;t--;}
         b[t]=a;   //4 jinzhi

//       for(int i=0;i<k;i++)
//       {
//           cout<<b[i]<<endl;
//       }

         // Count number of As Cs Gs Ts in a kmer word
         int count_ACGT[4]={0, 0, 0, 0};
         for (int mm = 0; mm < k; mm++)
         {
             count_ACGT[b[mm]]++;
         }
//         cout << "count " << count_ACGT[0] << count_ACGT[1] << count_ACGT[2] << count_ACGT[3] << endl;


         // Search specific kmer word probability
         long double pwREAL = TransToReal(Prior[DataNum][count_ACGT[0]][count_ACGT[1]][count_ACGT[2]])+TransToReal(Prior[DataNum][count_ACGT[3]][count_ACGT[2]][count_ACGT[1]]);

      
         // Output count-pw files              
         fout_CountPw << kk << "," << HashTableS[DataNum][kk] << "," << pwREAL << endl;
//         cout << kk << "," << HashTableS[DataNum][kk] << "," << pwREAL << endl;
         
         fout_CountPw.close();
              
      }


}

void D2C2compute(int k, unsigned long total[]) {
  int knew = k;
  for (int pp = 1; pp <= NS; pp++) {      
    c2Form[pp][pp] = 0; 
    c2starForm[pp][pp] = 0; 
    c2sheppForm[pp][pp] = 0; 
    c2QForm[pp][pp] = 0; 
    c2QstarForm[pp][pp] = 0; 
      
    for (int qq = pp + 1; qq <= NS; qq++) {                  
      SCIENTIFIC_NUMBER wordmean1, wordmean2;
      wordmean1.value = 0.0; wordmean1.factor = 0;
      wordmean2.value = 0.0; wordmean2.factor = 0;
      
      SCIENTIFIC_NUMBER tilde1, tilde2, tildeQ1, tildeQ2;
      tilde1.value = 0.0; tilde1.factor = 0;
      tilde2.value = 0.0; tilde2.factor = 0;
      tildeQ1.value = 0.0; tildeQ1.factor = 0;
      tildeQ2.value = 0.0; tildeQ2.factor = 0;

      
      SCIENTIFIC_NUMBER tildeprod, tildeprodQ;
      tildeprod.value = 0.0; tildeprod.factor = 0;
      tildeprodQ.value = 0.0; tildeprodQ.factor = 0;
      
      SCIENTIFIC_NUMBER vc;
      vc.value = 0.0; vc.factor = 0;

      SCIENTIFIC_NUMBER x2, y2, x2Q, y2Q;
      x2.value = 0.0; x2.factor = 0;
      y2.value = 0.0; y2.factor = 0;
      x2Q.value = 0.0; x2Q.factor = 0;
      y2Q.value = 0.0; y2Q.factor = 0;


      SCIENTIFIC_NUMBER x2tail, y2tail, x2Qtail, y2Qtail;
      x2tail.value = 0.0; x2tail.factor = 0;
      y2tail.value = 0.0; y2tail.factor = 0;
      x2Qtail.value = 0.0; x2Qtail.factor = 0;
      y2Qtail.value = 0.0; y2Qtail.factor = 0;


      SCIENTIFIC_NUMBER x2s, y2s;
      x2s.value = 0.0; x2s.factor = 0;
      y2s.value = 0.0; y2s.factor = 0;

      SCIENTIFIC_NUMBER d2, d2star, d2shepp, d2Q, d2Qstar;
      d2.value = 0.0; d2.factor = 0;
      d2star.value = 0.0; d2star.factor = 0;
      d2shepp.value = 0.0; d2shepp.factor = 0;
      d2Q.value = 0.0; d2Q.factor = 0;
      d2Qstar.value = 0.0; d2Qstar.factor = 0;

      SCIENTIFIC_NUMBER c2_deno, c2star_deno, c2shepp_deno, c2Q_deno, c2Qstar_deno;
      c2_deno.value = 0.0; c2_deno.factor = 0;
      c2star_deno.value = 0.0; c2star_deno.factor = 0;
      c2shepp_deno.value = 0.0; c2shepp_deno.factor = 0;
      c2Q_deno.value = 0.0; c2Q_deno.factor = 0;
      c2Qstar_deno.value = 0.0; c2Qstar_deno.factor = 0;

      SCIENTIFIC_NUMBER c2, c2star, c2shepp, c2Q, c2Qstar;
      c2.value = 0.0; c2.factor = 0;
      c2star.value = 0.0; c2star.factor = 0;
      c2shepp.value = 0.0; c2shepp.factor = 0;
      c2Q.value = 0.0; c2Q.factor = 0;
      c2Qstar.value = 0.0; c2Qstar.factor = 0;


      // Compute Pw for each Kmer from Prior[][][][]
      for(unsigned long kk = 0; kk < pow(4,knew); kk++) {
	int a; int b[30]; memset(b, 0, 30);
	a = kk; 
	int t=knew-1;
	while(a>=3) {b[t]=a%4;a/=4;t--;}
	b[t]=a;   //4 jinzhi

	//count A T C Gs in kmer word
	int count[4]={0, 0, 0, 0};
	for (int mm=0; mm<knew; mm++) {
	  count[b[mm]]++;
	}
	
	//search specific kmer word probability
	SCIENTIFIC_NUMBER pw1, pw2;
	
	// Pw_seq and Pw_inv_seq 
	pw1 = SciAddition(Prior[pp][count[0]][count[1]][count[2]],Prior[pp][count[3]][count[2]][count[1]]);
	pw2 = SciAddition(Prior[qq][count[0]][count[1]][count[2]],Prior[qq][count[3]][count[2]][count[1]]);
              

	wordmean1 = SciMultiple(pw1,total[pp]);
	wordmean2 = SciMultiple(pw2,total[qq]);

	SCIENTIFIC_NUMBER na_wordmean1, na_wordmean2;
	na_wordmean1.value = - wordmean1.value; na_wordmean1.factor = wordmean1.factor;
	na_wordmean2.value = - wordmean2.value; na_wordmean2.factor = wordmean2.factor;

	SCIENTIFIC_NUMBER Kmercount1, Kmercount2, KmercountQ1, KmercountQ2;

	Kmercount1 = TransToScientific(HashTableS[pp][kk]);
	Kmercount2 = TransToScientific(HashTableS[qq][kk]);
	KmercountQ1 = TransToScientific(HashTableSQ[pp][kk]);
	KmercountQ2 = TransToScientific(HashTableSQ[qq][kk]);




	// Average Word Probability as defined in qCluster [CLS14]
	SCIENTIFIC_NUMBER awp1, awp2, avg1, avg2;
	awp1 = TransToScientific(computeAWP(pp, kk));
	awp2 = TransToScientific(computeAWP(qq, kk));
	avg1 = SciMultiple(awp1, na_wordmean1);
	avg2 = SciMultiple(awp1, na_wordmean2);

	tilde1 = SciAddition(Kmercount1, na_wordmean1);
	tilde2 = SciAddition(Kmercount2, na_wordmean2);
	tildeQ1 = SciAddition(KmercountQ1, avg1);
	tildeQ2 = SciAddition(KmercountQ2, avg2);

	tildeprod = SciMultiple(tilde1, tilde2);
	tildeprodQ = SciMultiple(tildeQ1, tildeQ2);
	vc = SciPow(SciMultiple(wordmean1, wordmean2), 0.5);

	x2 = SciAddition(SciPow(Kmercount1,2), x2);
	y2 = SciAddition(SciPow(Kmercount2,2), y2);

	x2Q = SciAddition(SciPow(KmercountQ1,2), x2Q);
	y2Q = SciAddition(SciPow(KmercountQ2,2), y2Q);

	SCIENTIFIC_NUMBER inv_wordmean1, inv_wordmean2;
	inv_wordmean1.value = 1/wordmean1.value; inv_wordmean1.factor = - wordmean1.factor;
	inv_wordmean2.value = 1/wordmean2.value; inv_wordmean2.factor = - wordmean2.factor;
	x2tail = SciAddition(SciMultiple(SciPow(tilde1,2), inv_wordmean1), x2tail);
	y2tail = SciAddition(SciMultiple(SciPow(tilde2,2), inv_wordmean2), y2tail);
	// 
	x2Qtail = SciAddition(SciMultiple(SciPow(tildeQ1,2), inv_wordmean1), x2Qtail);
	y2Qtail = SciAddition(SciMultiple(SciPow(tildeQ2,2), inv_wordmean2), y2Qtail);


	SCIENTIFIC_NUMBER inv_temp1, temp1;
	temp1 = SciPow(SciAddition(SciPow(tilde1,2), SciPow(tilde2,2)), 0.5);
	inv_temp1.value = 1/temp1.value; inv_temp1.factor = - temp1.factor;
	x2s = SciAddition(SciMultiple(SciPow(tilde1,2),inv_temp1), x2s);
	y2s = SciAddition(SciMultiple(SciPow(tilde2,2),inv_temp1), y2s);

	d2 = SciAddition(SciMultiple(Kmercount1, Kmercount2), d2);
	d2Q = SciAddition(SciMultiple(KmercountQ1, KmercountQ2), d2Q);

	
	if(tildeprod.value != 0) {
	  SCIENTIFIC_NUMBER inv_temp2, temp2;
	  temp2 = SciPow(SciAddition(SciPow(tilde1,2), SciPow(tilde2,2)),0.5);
	  inv_temp2.value = 1/temp2.value; inv_temp2.factor = -temp2.factor;
	  d2shepp = SciAddition(SciMultiple(tildeprod, inv_temp2), d2shepp);
	  if(vc.value != 0) {
	    SCIENTIFIC_NUMBER inv_vc;
	    inv_vc.value = 1/vc.value; inv_vc.factor = - vc.factor;
	    d2star = SciAddition(SciMultiple(tildeprod, inv_vc), d2star);
	  }
	}

	if(tildeprodQ.value != 0) {
	  if(vc.value != 0) {
	    SCIENTIFIC_NUMBER inv_vc;
	    inv_vc.value = 1/vc.value; inv_vc.factor = - vc.factor;
	    d2Qstar = SciAddition(SciMultiple(tildeprodQ, inv_vc), d2Qstar);
	  }
	}

      }




      SCIENTIFIC_NUMBER inv_c2, inv_c2star, inv_c2shepp, inv_c2Q, inv_c2Qstar;
      c2_deno = SciMultiple(SciPow(x2, 0.5), SciPow(y2, 0.5));
      inv_c2.value = 1/c2_deno.value; inv_c2.factor = - c2_deno.factor;
      c2 = SciMultiple(d2, inv_c2);

      c2Q_deno = SciMultiple(SciPow(x2Q, 0.5), SciPow(y2Q, 0.5));
      inv_c2Q.value = 1/c2Q_deno.value; inv_c2Q.factor = - c2Q_deno.factor;
      c2Q = SciMultiple(d2Q, inv_c2Q);

            
      c2star_deno = SciMultiple(SciPow(x2tail, 0.5), SciPow(y2tail, 0.5));
      inv_c2star.value = 1/c2star_deno.value; inv_c2star.factor = - c2star_deno.factor;
      c2star = SciMultiple(d2star,inv_c2star);

      c2Qstar_deno = SciMultiple(SciPow(x2Qtail, 0.5), SciPow(y2Qtail, 0.5));
      inv_c2Qstar.value = 1/c2Qstar_deno.value; inv_c2Qstar.factor = - c2Qstar_deno.factor;
      c2Qstar = SciMultiple(d2Qstar, inv_c2Qstar);


      c2shepp_deno = SciMultiple(SciPow(x2s, 0.5), SciPow(y2s, 0.5));
      inv_c2shepp.value = 1/c2shepp_deno.value; inv_c2shepp.factor = - c2shepp_deno.factor;
      c2shepp = SciMultiple(d2shepp,inv_c2shepp);


      double d2Real, d2starReal, d2sheppReal, d2QReal, d2QstarReal, 
	c2Real, c2starReal, c2sheppReal, c2QReal, c2QstarReal;
      d2Real = TransToReal(d2);
      d2starReal = TransToReal(d2star);
      d2sheppReal = TransToReal(d2shepp);
      d2QReal = TransToReal(d2Q);
      d2QstarReal = TransToReal(d2Qstar);
      c2Real = TransToReal(c2);
      c2starReal = TransToReal(c2star);
      c2sheppReal = TransToReal(c2shepp);
      c2QReal = TransToReal(c2Q);
      c2QstarReal = TransToReal(c2Qstar);


      c2Form[pp][pp] = 0; c2Form[qq][pp] = c2Form[pp][qq] = 1.0 - c2Real;
      c2starForm[pp][pp] = 0; c2starForm[qq][pp] = c2starForm[pp][qq] = 1.0 - c2starReal;
      c2sheppForm[pp][pp] = 0; c2sheppForm[qq][pp] = c2sheppForm[pp][pp] = 1.0 - c2sheppReal;
      
      c2QForm[pp][pp] = 0; c2QForm[qq][pp] = c2QForm[pp][qq] = 1.0 - c2QReal;
      c2QstarForm[pp][pp] = 0; c2QstarForm[qq][pp] = c2QstarForm[pp][qq] = 1.0 - c2QstarReal;
      
      d2Form[pp][pp] = 0; d2Form[qq][pp] = d2Form[pp][qq] = d2Real;
      d2starForm[pp][pp] = 0; d2starForm[qq][pp] = d2starForm[pp][qq] = d2starReal;

      d2QForm[pp][pp] = 0; d2QForm[qq][pp] = d2QForm[pp][qq] = d2QReal;
      d2QstarForm[pp][pp] = 0; d2QstarForm[qq][pp] = d2QstarForm[pp][qq] = d2QstarReal;

      
    }
    
  }

}




void PrintD2C2(int d, char* argv_DataNumK)
{

     char dstr[5]; sprintf(dstr, "%d", d); 
       
     // output csv file
     char outputfile_c2[1000];
     strcpy(outputfile_c2,"k");
     strcat(outputfile_c2,argv_DataNumK);
     strcat(outputfile_c2,"-d");
     strcat(outputfile_c2,dstr);
     strcat(outputfile_c2,"-c2Form");
//     cout << outputfile_c2 << endl;
     ofstream fout_c2(outputfile_c2);

     for(int t1 = 1; t1 <= NS; t1++)
     {
          for(int t2 = 1; t2 <= NS; t2++)
          {
	    if(t2 <= t1) { 
	      fout_c2 << c2Form[t1][t2] << "," ; 	       

	    }
              
	    else if(t2 >t1 && t2 != NS ) {
	      fout_c2 << ","; 

	    }
               else {
		 fout_c2 << "\n"; 

	       }  

           }
     }
       
     fout_c2.close();




      char outputfile_c2star[1000];
      strcpy(outputfile_c2star,"k");
      strcat(outputfile_c2star,argv_DataNumK);
      strcat(outputfile_c2star,"-d");
      strcat(outputfile_c2star,dstr);
      strcat(outputfile_c2star,"-c2starForm");
      ofstream fout_c2star(outputfile_c2star);
           
      for(int t1 = 1; t1 <= NS; t1++)
      {
           for(int t2 = 1; t2 <= NS; t2++)
           {
                if(t2 <= t1) { fout_c2star << c2starForm[t1][t2] << "," ; }
                if(t2 >t1 && t2 != NS ) { fout_c2star << ","; }
                if(t2==NS){fout_c2star << "\n"; }                     
           }
      }
       
      fout_c2star.close();



     char outputfile_c2shepp[1000];
     strcpy(outputfile_c2shepp,"k");
     strcat(outputfile_c2shepp,argv_DataNumK);
     strcat(outputfile_c2shepp,"-d");
     strcat(outputfile_c2shepp,dstr);
     strcat(outputfile_c2shepp,"-c2sheppForm");
     ofstream fout_c2shepp(outputfile_c2shepp);


     for(int t1 = 1; t1 <= NS; t1++)
     {
          for(int t2 = 1; t2 <= NS; t2++)
          {
                    
               if(t2 <= t1) {  fout_c2shepp <<c2sheppForm[t1][t2] << "," ; }
               if(t2 >t1 && t2 != NS ) { fout_c2shepp << ","; }
               if(t2==NS){fout_c2shepp << "\n"; }                     

          }
     }
       
     fout_c2shepp.close();


}





// Extract the raw data(.fasta/fastq) file name and its direction. 
int AllocateFile(char* argv_DataFileName){

        char seq[500][999];
        int seqlength=999; 
        char delims[] = " ";
        char *dir_name = NULL;
      
 
        ifstream infile(argv_DataFileName); 
        int lineNum = 1;

        while(infile.getline(seq[lineNum],seqlength))
        {	

            dir_name = strtok( seq[lineNum], delims );

            int i = 0;
            while( dir_name != NULL ){
                
                if(i == 0){
                    
                    fileDirect[lineNum] = dir_name; 
                         
                    dir_name = strtok( NULL, delims );
                    i ++;
                    continue; 
                                      
                }
                
                if(i == 1){
                        
                    fileName[lineNum] = dir_name; 
                    break;    
                }
            }

            
            lineNum ++;
            

   
       }
    
       

       return (lineNum-1);


}
    
    

// Data Pre-process: Transform .fasta(/.fastq) to read sequence only data.
void dataPreProcess(int DatasetNum, char *fastType){


    char outputfile[1000];
    strcpy(outputfile,fileName[DatasetNum]);
    strcat(outputfile,"_combine");
    ofstream fout(outputfile); 
 

    char seq[99999];
    int seqlength=99999;
    ifstream infile(fileDirect[DatasetNum]);
    

    if(fastType == "Q" || fastType == "q"){
         
         cout << seq << endl;
        
         int flagA = 0;
         while(infile.getline(seq,seqlength))
         {
             
             
    	       if(seq[0]=='>')
    	       {
    	             if(flagA == 1){fout << endl;}
                     continue;
    	       }
    	       else
    	       {
                     flagA = 1;
                     fout << seq;
               }
         }
    }
    else{
         
         int flagQ = 0;
        
         while(infile.getline(seq,seqlength))
         {
             
    	       if(seq[0]=='A' || seq[0]=='a' || seq[0]=='C' || seq[0]=='c' || seq[0]=='G' || seq[0]=='g' || seq[0]=='T' || seq[0]=='t' || seq[0]=='N' || seq[0]=='n')
    	       {
                     flagQ = 1;
    	             fout<<seq;
    	       }
    	       else
    	       {
                     if(flagQ == 1){fout << endl;}
                     flagQ = 0;
                     continue;
               }
         }  


   }

}

void printOutMatrix(double m[200][200]) {
  for (int i = 1; i <= NS; ++i) {
    for (int j = 1; j <= NS; ++j) {
      std::cout << m[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

void printAllMatrices() {
  std::cout << "\n-- D2 --" << std::endl;
  printOutMatrix(d2Form);
  std::cout << "\n-- D2star --" << std::endl;
  printOutMatrix(d2starForm);
  std::cout << "\n-- d2 --" << std::endl;
  printOutMatrix(c2Form);
  std::cout << "\n-- d2star --" << std::endl;
  printOutMatrix(c2starForm);    
  std::cout << "\n-- D2Q --" << std::endl;
  printOutMatrix(d2QForm);
  std::cout << "-- D2Qstar --" << std::endl;
  printOutMatrix(d2QstarForm);  
  std::cout << "\n-- d2Q --" << std::endl;
  printOutMatrix(c2QForm);
  std::cout << "-- d2Qstar --" << std::endl;
  printOutMatrix(c2QstarForm);

}

void printHashtables() {
  for (int d = 1; d <= NS; ++d) {
    std::cout << "H(" << d <<  "):  ";
    for( hash_map<unsigned long,unsigned long>::iterator i = HashTableS[d].begin(); i!= HashTableS[d].end(); i++) {
      std::cout << i->second << " ";
    }
    std::cout << std::endl;
  }
  for (int d = 1; d <= NS; ++d) {
    std::cout << "Hq(" << d << "): ";
    for( hash_map<unsigned long,double>::iterator i = HashTableSQ[d].begin(); i!= HashTableSQ[d].end(); i++) {
      std::cout << i->second << " ";
    }
    std::cout << std::endl;
  }
}

void printAWPs() {
  for (int d = 1; d <= NS; ++d) {
    std::cout << "AWP(1): [ ";
    unsigned long K = pow(4,k);
    for (unsigned long index = 0; index < K; ++index) {
      std::cout << computeAWP(d,index) << " ";
    }
    std::cout << "]" << std::endl;
  }
}


void printFilesInfo() {
  for (int d = 1; d <= NS; ++d) {
    std::cout << d << " " << fileName[d] << "   " << fileDirect[d] << endl;
  }
}


void printPrior() {
  std::cout << "\n-- Pinting Prior --" << endl;
  for (int d = 1; d <= NS; ++d) {
    std::cout << d << ". ";
    for (int a = 0; a < k+1; ++a) {
      for (int c = 0; c < k+1; ++c) {
	for (int g = 0; g < k+1; ++g) {
	  cout << TransToReal(Prior[d][a][c][g]) << " ";
	}
      }
    }
    std::cout << std::endl;
  }  
}

void writeMatrixToFile(double m[200][200], const std::string& fileName) {
  std::ofstream fout(fileName.c_str());
  for (int i = 1; i <= NS; ++i) {
    for (int j = 1; j <= NS; ++j) {
      fout << m[i][j] << "\t";
    }
    fout << std::endl;
  }
  fout.close();
}

void writeMatricesToFile() {
  writeMatrixToFile(c2Form, "d2.out");
  writeMatrixToFile(c2starForm, "d2star.out");

  writeMatrixToFile(c2QForm, "d2Q.out");
  writeMatrixToFile(c2QstarForm, "d2Qstar.out");

  writeMatrixToFile(d2Form, "D2.out");
  writeMatrixToFile(d2starForm, "D2star.out");

  writeMatrixToFile(d2QForm, "D2Q.out");
  writeMatrixToFile(d2QstarForm, "D2Qstar.out");
}

void printParameters() {
  std::cout << "\n------------------------ PARAMETERS ------------------------" << std::endl;
  std::cout << "*  k                  "  << k << std::endl;
  std::cout << "*  NS                 "  << NS << std::endl;
  std::cout << "*  M                  "  << M<< std::endl;
  std::cout << "*  type               "  << type << std::endl;
  std::cout << "------------------------------------------------------------\n" << std::endl;
}

void printUsage() {
  std::cout << endl;
  std::cout << "USAGE:" << std::endl;
  std::cout << "   c2q [k] [list] [type] [M]" << std::endl;
  std::cout << std::endl;
  std::cout << "      k        k-mer length" << std::endl;
  std::cout << "      list     a file containint list of input sequences file in the form [fast file] [seq name]" << std::endl;
  std::cout << "      type     either Q for quality value based ouput or A for non quality output" << std::endl;
  std::cout << "      M        length of the reads (needed to compute averages)" << std::endl;    
  std::cout << endl;
}

//KmerCount.out [sample data file] [k]
int main(int argc, char *argv[]) {

if (argc < 5) {
  printUsage();
  return 1;
 }

// Extract the User-specific parameters.
 k = atoi(argv[1]);
 NS = AllocateFile(argv[2]);    
 char *dataType = argv[3];   
 type = *dataType;
 M = atoi(argv[4]);
 printParameters();

 HashTable.clear();
 HashPw.clear();
 unsigned long total[NS+1];     //*****change parameters

//    cout << fileName[1] << fileName[2] << fileName[3] << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Data Pre-process: Eliminate the annotation lines in raw dataset .fastq/.fasta. Only leave the read sequences
    cout << "0. Data Pre-process: Eliminate the annotation lines in raw dataset .fastq/.fasta. Only leave the read sequences..." << endl;
    
    for (int DataNum = 1; DataNum <= NS; DataNum++)     //////
    {
          
        cout << "Data Pre-process Dataset No." << DataNum << "." << endl;
        
	        dataPreProcess(DataNum, dataType);
 
    }





    // Count frequency for each Kmer word and compute Pw for each Kmer word in each dataset
    cout << "1. Count frequency for each Kmer word and compute Pw for each Kmer word in each dataset..." << endl;
    
    for (int DataNum = 1; DataNum <= NS; DataNum++)     //////
    {
          
        cout << "Scanning Dataset No." << DataNum << "." << endl;
        
         // Count kmer in rth dataset
         total[DataNum] = DataKmerCount(DataNum, fileDirect[DataNum], k);

         // Compute Probability of each kmer word(Pw)
         ComputePw(DataNum, k, argv[1], fileName[DataNum], d);
 
    }




    // Output key, hashtable and pw for each dataset
    
    cout << "2. Print out key, hashtable and pw for each dataset..." << endl;
 
    for (int DataNum = 1; DataNum <= NS; DataNum++)  
    {
        
         cout << "Printing out Kmer count and its Pw for Dataset No." << DataNum << "." << endl;

         // Print kmer count hashtable and its corresponding Pw
         PrintKmerCountPw(k, DataNum, argv[1],fileName[DataNum], d);

    }



    // Compute D2C2 for any pair of datasets from their Kmer-Pw information 
    
    cout << "3. Compute D2C2 for any pair of datasets from their Kmer-Pw information..." << endl;

//    cout << "k0" << k << endl;
    
    D2C2compute(k, total);


    
    // Output D2C2 result matrix: c2Form, c2starForm, c2sheppForm
    
    cout << "4. Print D2C2 result matrix: c2Form, c2starForm, c2sheppForm..." << endl;
    writeMatricesToFile();


    //    PrintD2C2(d, argv[1]);

    //printPrior();
    // printFilesInfo();
    //printHashtables();
    //printAllMatrices();
    
    // std::cout << "\n-- d2star --" << std::endl;
    // printOutMatrix(c2starForm);    
    // std::cout << "-- d2Qstar --" << std::endl;
    // printOutMatrix(c2QstarForm);


    return 0;
}






/////////////////////////////////////END/////////////////////////////////////////////////////

