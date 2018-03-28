#pragma once 
#include "network.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <string>

#define mu	1.0                                                     //�A�Ȳv
#define SF   1.0                                                    //Spreading factor�A�T�w
#define Pn  (double) 4.11*pow(10.0,(-14))							//�����T(mW)�A�`�N���A�䤤	Pn=kT k��Boltzmann�`��	T������ū�
#define Keq 2.75*pow(10.0,(15))										//Gain�����̷|�Ψ�
#define Alpha (float) 1.0											//����]�l�A�T�w
#define W  10.0														//�W�e , �p�G�o���F���`�N�����TPn�]�n��
#define Ps 200.0													//��a�x�R�A�\�vW
#define bsmaxpower 40.0												//��a�x�үണ�ѥ\�v�`�M�W��W
#define k    1.0                                                    //�ϥΪ̳��ǿ�q�һݤ�I���O��(dollars/bit)      //#define k  0.2
#define c    1.0                                                    //����q�q�O(dollars/joule)          //#define c	13.0/1000
#define pi   5.0                                                    //���Q��n���h�֩ҳ]���ܼ�
#define Big  (double) 20.0											//Big M
#define dmax 1500.0													//�̤j�Z��
#define beta  1.0													//�ؼЦ�	

#define ga 6.0														//1�����ä��q�A3~10 ���D���ä��q(�V�j�N����V��)

#define par (int) 4													//bilinear term���q(���V�h�V����)##
#define nsum 4.0

//�`�ΰѼƭק�
#define PUB  600.0													//Power�̤j(mW)
#define SINR   100.0												//SINR�̤j
#define Rate 0.8													//datarate



using namespace std; 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
 

void NETWORK :: AddNode_T(NODE *node){                                           //�[�J�s��TP
	if(Node_T==NULL){                                                            //�p�G�쥻���O�Ū��A�N�O�Ĥ@��
		Node_T = new LIST<NODE*>(node);                                          //new�B��l�|�t�m�@��data type�һݭn���Ŷ��A�öǦ^�ӪŶ�����}�F�ϥ�new�B��l�ʺA�t�m���Ŷ��A�b��ӵ{�������e�ä��|�۰��k�ٵ��O����A�����ϥ�delete�N�o�ӪŶ��ٵ��O��
		Node_T->GetValue()->SetNo(N++);                                          //��ȫ�A�Ǹ�++
		//cout<<"n="<< Node_T->GetValue()-> GetNo() <<endl;
		 
	}  else{                                                                     //�p�G���O�Ū�
		Node_T->Tail()->SetNext(new LIST<NODE*>(Node_T->Tail(),node,NULL));      //�q���ڥ[�J�A�̷�(Node_T->Tail(),node,NULL)������
		Node_T->Tail()->GetValue()->SetNo(N++);                                  //��ȫ�A�Ǹ�++
	} 
}

void NETWORK :: AddNode_B(NODE *node){                                           //�y�k�PTP�ۦP 
	if(Node_B==NULL){
		Node_B = new LIST<NODE*>(node);
		Node_B->GetValue()->SetNo(M++);
		 
	} else{
		Node_B->Tail()->SetNext(new LIST<NODE*>(Node_B->Tail(),node,NULL));
		Node_B->Tail()->GetValue()->SetNo(M++);
	} 
} 




//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
void NETWORK :: AddUBLB(){                         //�W�B�U���禡

int M,N,n,n2,n3,n4;                                //��a�x�ƥ�M�A�ϥΪ̼ƥ�N

	LIST<NODE*> *node,*node2;                      //i,j
	LIST<NODE*> *node3,*node4;                     //i',j'
	

	M = this->GetNumberOfNode_B();                 //���oBS�ƶq	 
	N = this->GetNumberOfNode_T();                 //���oTP�ƶq

		//cout<<"UBLBHERE"<<endl;
		//system("pause");


//----------------------------------------------------------------------��l��-------------------------------------------------------------------------------------------------------//


	this->BinVar_LB = new int[ M + M*N + N + 2*M*N*par ];			 //Binary Variable�AM��yi�AM*N��xi,j�AM*N*par���fv,�fk i',j',np�AN�Ozj
	this->BinVar_UB = new int[ M + M*N + N + 2*M*N*par ];
	for(int i=0;	  i< M + M*N + N + 2*M*N*par  ; i++){            //�ΰj���l�ƩҦ��O�����}
		BinVar_LB[i] = 0;
		BinVar_UB[i] = 0;
	}

//------------------------------------------------------------------------------------------------

	this->Var_LB = new float[ 1 + 2*M*N + 3*N + M + 2*M*N*par ] ;         //Power,SINR,H=(1+SINR),�u��Downlink�A�ЫذʺA�O����Ŷ�(����)
	this->Var_UB = new float[ 1 + 2*M*N + 3*N + M + 2*M*N*par ] ;

	for(int i=0;   i< 1 + 2*M*N + 3*N + M + 2*M*N*par ; i++){             //�ΰj���l�ƩҦ��O�����}
		Var_LB[i] = 0; 
		Var_UB[i] = 0;
	}

//-------------------------------------------------------------------------------------------------

	this->BilVar_LB = new float[  2*N*M*N + 2*N*M*N*par ]; //Biliear Tearm�A2*N*M*N��K & V�A�F2*N*M*N*par���GK & �GV�Apar�����O��n�q
	this->BilVar_UB = new float[  2*N*M*N + 2*N*M*N*par ];

	for(int i=0;      i< 2*N*M*N + 2*N*M*N*par ; i++){     //�ΰj���l�ƩҦ��O�����}
		 BilVar_LB[i] = 0; 
		 BilVar_UB[i] = 0; 
	}
 
//------------------------------------------------------------�ܼƤW�U��------------------------------------------------//
 //u
 
     cout<<"Connection = 0.8"<<endl;
     cout<<"Rate = 0.35M"<<endl;
     cout<<"Ps = 200W "<<endl;
    

	//Var_UB[ 1+ M + 4*M*N ] =1/(20*M*N + (Ps*M));				//u ���ƭp�⪺����

	 Var_UB[ 2*M*N + 3*N + M  ] = 2.0;							//u��Upper Bound
    //cout<<"u="<<Var_UB[ 2*M*N + 3*N + M  ]<<endl;

     Var_LB[ 2*M*N + 3*N + M  ] = 0;							//u��Lower Bound


//------------------------------------------------------------------------Binary Term && Variable Term---------------------------------------------//

node = this-> GetNode_B();                                                           //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                                   //i

	node2  = this-> GetNode_T();                                                     //��node����}����TP
	while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                                             //j
				
//Varible �ŧi
						
		
		        Var_UB[ 	    n*N + n2  ] = 250 * Var_UB[ 2*M*N + 3*N + M  ];                           //Power
				Var_UB[               M*N + n2  ] = SINR  * Var_UB[ 2*M*N + 3*N + M  ];	                  //SINR*u 
				Var_UB[           M*N + N + n2  ] = log(1+Var_UB[ M*N + n2 ]);                            //H = (1+SINR)		
                Var_UB[   M*N + 2*N + n*N + n2  ] = Var_UB[ 2*M*N + 3*N + M ];                            //x*u = u
                Var_UB[       2*M*N + 2*N + n   ] = Var_UB[ 2*M*N + 3*N + M ];                            //y*u = u
				Var_UB[   2*M*N + 2*N + M + n2  ] = Var_UB[ 2*M*N + 3*N + M ];                            //Zj*u = u
			//	Var_UB[       2*M*N + 3*N + M   ] = Var_UB[ 2*M*N + 3*N + M ];                            //u





//Binary �ŧi		
		        BinVar_UB [           n   ] = 1.0;                                            //yi
				BinVar_UB [ M + n*N + n2  ] = 1.0;                                            //xij
				BinVar_UB [ M + M*N + n2  ] = 1.0;                                            //Zj 			  



			node2  = node2  ->GetNext();	
		} 
		 node  = node  ->GetNext();	
	}



//------------------------------------------------------------------------Bilinear Term---------------------------------------------//

//( V & K )*U

	 		 
		node = this -> GetNode_B();																				//��node����}����BS
	while ( node != NULL )
	{
		n = node -> GetValue() -> GetNo();																		//i
		node2 = this -> GetNode_T();																			//��node����}����TP
		
		while ( node2 != NULL )
		{
			n2 = node2 -> GetValue() -> GetNo();																//j
	 		node3 = this -> GetNode_B(); 
			
			while ( node3 != NULL )
			{
				n3 = node3 -> GetValue() -> GetNo();															//i'
				node4 = this -> GetNode_T();
				
				while ( node4 != NULL )
				{
					n4 = node4 -> GetValue() -> GetNo();														//j'


			                                   BilVar_UB[     n2*M*N + n*N + n4 ] = 100 * SINR * Var_UB[ 2*M*N + 3*N + M ] ;                           //V�����쪺�O�����m => Powerij'*SINRj*U
					                       
											   BilVar_UB[     N*M*N + n2*M*N + n3*N + n4 ] = 100 * SINR * Var_UB[ 2*M*N + 3*N + M ] ;                  //K�����쪺�O�����m => Poweri'j'*SINRj*U	                               
				 
				 node4 = node4 -> GetNext();	
		      } 
			  node3 = node3 -> GetNext();		
		   } 
			node2 = node2 -> GetNext();	
		} 
		node  = node -> GetNext();	
	} 
//��V & ��K


	node = this -> GetNode_B();																					//�P�W
	while ( node != NULL )
	{
		n = node -> GetValue() -> GetNo();																		//i
		node2 = this -> GetNode_T();
		
		while ( node2 != NULL )
		{
			n2 = node2 -> GetValue() -> GetNo();																//j
	 		node3 = this -> GetNode_B(); 
			
			while ( node3 != NULL )
			{
				n3 = node3 -> GetValue() -> GetNo();															//i'
				node4 = this -> GetNode_T();
				
				while ( node4 != NULL )
				{
					n4 = node4 -> GetValue() -> GetNo();														//j'
			        
					        
					
	              for(int j=0; j<par; j++){
					                           	 BilVar_UB[                 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ] = PUB * SINR * Var_UB[ 2*M*N + 3*N + M ] ;           //��V = Powerij'*SINRj*U��upper bound

												 BilVar_UB[     2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ] = PUB * SINR * Var_UB[ 2*M*N + 3*N + M ] ;         //��K = Poweri'j'*SINRj*U��upper bound	 
  
				                         
				                              
			    
				  }
				 node4 = node4 -> GetNext();	
		      } 
			  node3 = node3 -> GetNext();		
		   } 
			node2 = node2 -> GetNext();	
		} 
		node  = node -> GetNext();	
	} 
	

//-----------------------------------------------------------Partition---------------------------------------------------

//�fv && �fv*u
node = this-> GetNode_B();                                                            //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                                    //i

	node4  = this-> GetNode_T();                                                      //��node����}����TP
	    while(node4 != NULL){
		n4 = node4->GetValue()->GetNo(); 	                                          //j'	                                          

             for(int j=0; j<par; j++)
			 {
                                 
				   BinVar_UB[       M + M*N + N + n*N*par + n4*par + j  ] = 1 ;                                         //�N�b�ۦP��n��n4�U���q����j�q�]���fv
				   Var_UB[  1 + 2*M*N + 3*N + M + n*N*par + n4*par + j  ] = Var_UB[ 2*M*N + 3*N + M ];                  //�fv*u = u			   
			 } 
		  
		node4 = node4 -> GetNext();	                                                  //���V�U�@��TP
	   }  
	node = node -> GetNext();	                                                      //���V�U�@��BS
   } 

//�fk  && �fk*u
 node3 = this-> GetNode_B();                                                           //��node����}����BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();                                                   //i'

	node4 = this-> GetNode_T();                                                        //��node����}����TP
	    while(node4 != NULL){
		n4 = node4->GetValue()->GetNo(); 	                                           //j'	                                           

             for(int j=0; j<par; j++)
			 {
                                 		   
				   BinVar_UB[       M + M*N + N + M*N*par + n3*N*par + n4*par + j  ] = 1 ;                               //�N�b�ۦP��n3��n4�U���q����j�q�]���fk
				   Var_UB[  1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j  ] = Var_UB[ 2*M*N + 3*N + M ];        //�fk*u = u
			 } 
		        
     	node4 = node4 -> GetNext();	                                                   //���V�U�@��TP
	   }  
	node3 = node3 -> GetNext();	                                                       //���V�U�@��BS
   }

    return;
	system("pause");

	}

//----------------------------Add Variables ��m---------------------------------------------------------//

void NETWORK ::func(){ //NETWORK *network

	//cout<<"Vars Location"<<endl;


    int M,N;

	M = this->GetNumberOfNode_B();	 
	N = this->GetNumberOfNode_T();
	
	ILOSTLBEGIN                                                          //cplex���зǨ禡�w
 
	
	IloEnv env;                                                          //�ظmenv����
	IloModel model(env);                                                 //env���Ҹ̫إ�model�ҫ�
	IloNumVarArray vars(env), Vars(env), binvars(env), bilvars(env);     //�ŧi�}�C
	IloNumArray Coes(env), vals(env);                                    //�ŧi�ƭ�
	
	

	LIST<NODE*> *node, *node2;                                          //i,j
	LIST<NODE*> *node3,*node4;                                          //i',j'
	LIST<NETWORK*> *network;											//�o�̭��ƫŧi�O���F�N��}�ɤJcplex

	int n,n2,n3,n4,i,j;	

 
  try{                                                                   //���~�����y�k

		//if(BNB_DDisplay){
		//	cout << "Adding vars."<<endl;

	
 char var_name[60];                                                      //�]�@�Ӥj�p��60���}�C,�s�ܼƦW�r
// char cst_name[60];                                                    //�]�@�Ӥj�p��60���}�C,�s�`�ƦW�r
		 

// Add Power_D
node  = this-> GetNode_B();                                              //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                       //i

    node2  = this-> GetNode_T();                                         //��node2����}����TP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo();                             //j
				
				char var_name[60];                                                                             //�]�@�Ӥj�p��60���}�C ##
				sprintf(var_name,"P_%d_%d",n+1,n2+1);                                                          //��P_D_%d_%d�g�Jvar_name�A�æL�X�ӰѼ�i,j
			 	vars.add(IloNumVar(env, Var_LB[ n*N + n2 ], Var_UB[ n*N + n2 ], var_name));                    //��c�[��cplex�Apower_D�n�b�W�U������

					  
			 node2 = node2 -> GetNext();	
		    } 
			 
 
		node  = node  -> GetNext();	
		}  


//Add SINR_D                                           
	
    node2  = this-> GetNode_T();                                           //�PPower_D�y�k�ۦP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo(); 
			
				char var_name[60];                                                                              
				sprintf(var_name,"S_%d" ,n2+1);                                                                
			 	 vars.add(IloNumVar(env, Var_LB[ M*N + n2 ], Var_UB[ M*N + n2 ], var_name));
				
				 node2 = node2 -> GetNext();	
		        }
	 
			  

//Add H_D = 1+SINR

    node2  = this-> GetNode_T();                                           //�PPower_D�y�k�ۦP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo(); 

			char var_name[60];
			sprintf(var_name,"H_%d",n2+1);
		 	vars.add(IloNumVar(env, Var_LB[ M*N + N + n2 ], Var_UB[ M*N + N + n2 ], var_name));
 
			     node2 = node2 -> GetNext();	
		        }
			 

//Add x*u

node  = this-> GetNode_B();																			//��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();																	//i

    node2  = this-> GetNode_T();																	//��node����}����TP
		while (node2 != NULL){  
	    n2 = node2->GetValue()->GetNo();															//j

			char var_name[60];
			sprintf(var_name,"x_%d_%d",n+1,n2+1);
		 	vars.add(IloNumVar(env, Var_LB[ M*N + 2*N + n*N + n2 ], Var_UB[ M*N + 2*N + n*N + n2 ], var_name));   //x*u���O�����m
			 
	
			     node2 = node2 -> GetNext();	
		        } 
			 
	 
			 node  = node  -> GetNext();	
		    }  

//Add y*u

node  = this-> GetNode_B();																			//��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();																	//i

			char var_name[60];
			sprintf(var_name,"y_%d",n+1 );
		 	vars.add(IloNumVar(env, Var_LB[ 2*M*N + 2*N + n ], Var_UB[ 2*M*N + 2*N + n ], var_name));                  //y*u���O�����m
 
			 node  = node  -> GetNext();	
		    }  
 
//Add Zj*u

node2  = this-> GetNode_T();																		//��node����}����BS
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();																//j

			char var_name[60];
			sprintf(var_name,"Zj_%d",n2+1 );
		 	vars.add(IloNumVar(env, Var_LB[ 2*M*N + 2*N + M + n2 ], Var_UB[ 2*M*N + 2*N + M + n2 ], var_name));        //Zj*u���O�����m
 
			 node2  = node2  -> GetNext();	
		    } 

	//Add u
	
	//char var_name[60];
	sprintf(var_name,"u");
	vars.add(IloNumVar(env, Var_LB[ 2*M*N + 3*N + M ], Var_UB[ 2*M*N + 3*N + M ], var_name));                          //u���O�����m


//Add �fv*u

node  = this-> GetNode_B();                                                                         //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();																	//i

    node4  = this-> GetNode_T();																	//��node����}����TP
		while (node4 != NULL){  
			n4 = node4->GetValue()->GetNo();														//j'

				for(int j=0 ; j<par ; j++)
				{
					char var_name[60];
					sprintf(var_name,"�fv_%d_%d_%d",n+1,n4+1,j+1);
		 			vars.add(IloNumVar(env, Var_LB[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ], Var_UB[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ], var_name));                       //�fv*u���O�����m
			    }
	
			  node4 = node4 -> GetNext();	
		     }  
			 
 
		node  = node  -> GetNext();	
	   }  

//Add �fk*u

node3  = this-> GetNode_B();                                                                          //��node����}����BS
	while(node3 != NULL){                                                                             
	n3 = node3->GetValue()->GetNo();                                                                  //i'

    node4  = this-> GetNode_T();                                                                      //��node����}����TP
		while (node4 != NULL){  
			n4 = node4->GetValue()->GetNo();                                                          //j'

				for(int j=0 ; j<par ; j++)
				{
					char var_name[60];
					sprintf(var_name,"�fk_%d_%d_%d",n3+1,n4+1,j+1);
		 			vars.add(IloNumVar(env, Var_LB[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ], Var_UB[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ], var_name));    //�fk*u���O�����m
			    }
	
			  node4 = node4 -> GetNext();	
		     }  
			 
 
	      node3  = node3  -> GetNext();	
		 }




//------------------------------------------------------Add Binary ��m----------------------------------------------------------------------//

//Add y
node  = this-> GetNode_B();                                                                           //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                                                    //i

	      		char var_name[60];                                                                    //�ŧi�ܼƦW�٪��O���� ##
				sprintf(var_name,"Bin_y_%d",n+1);                                                     //��Bin_y_%d�g�Jvar_name�A�æL�Xi,j
			 	binvars.add(IloNumVar(env, BinVar_LB[ n ], BinVar_UB[ n ],ILOINT , var_name));        //�Ny�[�Jcplex���Ҥ�
				//cout <<"y:BinVar_LB = "<<i<< endl;
		   	  node  = node  -> GetNext();	
		    } 
 
//Add x
node  = this-> GetNode_B();                                                                                       //��node����}����BS
	while(node != NULL){                                                                         
	n = node->GetValue()->GetNo();                                                                                //i
 
    node2  = this-> GetNode_T();                                                                                  //��node����}����TP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo();                                                                      //j

			char var_name[60];                                                                                    //�ŧi�ܼƦW�٪��O���� ##
			sprintf(var_name,"Bin_x_%d_%d",n+1,n2+1);                                                             //��Bin_x_%d_%d�g�Jvar_name�A�æL�Xi,j
		 	binvars.add(IloNumVar(env, BinVar_LB[ M + n*N + n2 ], BinVar_UB[ M + n*N + n2 ],ILOINT , var_name));  //�Nx�[�Jcplex���Ҥ�
 
			 
			   node2 = node2 -> GetNext();	
		      }
			 
		   node  = node  -> GetNext();
	      }

	
//Add Zj
node2  = this-> GetNode_T();                                                                                                //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                                                                        //j

	      		char var_name[60];                                                                                          //�ŧi�ܼƦW�٪��O���� ##
				sprintf(var_name,"Bin_Zj_%d",n2+1);                                                                         //��Bin_Zj_%d�g�Jvar_name�A�æL�Xj
			 	binvars.add(IloNumVar(env, BinVar_LB[ M + M*N + n2 ], BinVar_UB[ M + M*N + n2 ],ILOINT , var_name));        //�NZj�[�Jcplex���Ҥ�
	
		   	  node2  = node2  -> GetNext();	
		    } 
				
//Add�fv  	  
node  = this-> GetNode_B();                                                                                                                                        
	while(node != NULL){ 
	n = node->GetValue()->GetNo();                                                                                
 
    node4  = this-> GetNode_T();                                                                                   
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                                                       

						for(int j=0 ; j<par ; j++)
						{

							char var_name[60];
							sprintf(var_name,"Bin_�fv_%d_%d_%d ",n+1,n4+1,j+1);                                                                                                        //��Bin_�fv_%d_%d_%d�g�Jvar_name�A�æL�Xi,j,par
		 					binvars.add(IloNumVar(env, BinVar_LB[ M + M*N + N + n*N*par + n4*par + j ], BinVar_UB[ M + M*N + N + n*N*par + n4*par + j ],ILOINT , var_name));           //�N�fv�[�Jcplex���Ҥ�
						   
						}
 				  node4 = node4-> GetNext();
				 }   
	        
			node = node -> GetNext();
	       }

//Add�fk  	  
node3  = this-> GetNode_B();                                                                                                                                        
	while(node3 != NULL){ 
	n3 = node3->GetValue()->GetNo();                                                                                
 
    node4  = this-> GetNode_T();                                                                                   
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                                                       

						for(int j=0 ; j<par ; j++)
						{

							char var_name[60];
							sprintf(var_name,"Bin_�fk_%d_%d_%d ",n3+1,n4+1,j+1);                                                                                                                        //��Bin_�fk_%d_%d_%d�g�Jvar_name�A�æL�Xi,j,par
		 					binvars.add(IloNumVar(env, BinVar_LB[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ], BinVar_UB[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ],ILOINT , var_name));      //�N�fk�[�Jcplex���Ҥ�
						   
						}
 				  node4 = node4-> GetNext();
				 }   
	        
			node3 = node3 -> GetNext();
	       }

 
//-----------------------------------------------------Add Bilinear Term ��m------------------------------------------------------//
 
//Add V=Powerij'*SINRj*U		

	node2  = this-> GetNode_T();                                                                                                                                //��node����}����TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                                                                                                                        //j
			//cout<< "n2=" << node2->GetValue()->GetNo() <<endl;
 		   node = this-> GetNode_B();                                                                                                                           //��node����}����BS
			  while(node != NULL){
			  n = node->GetValue()->GetNo();                                                                                                                    //i
				//cout<< "n3=" << node3->GetValue()->GetNo() <<endl;
                 node4 = this-> GetNode_T();                                                                                                                    //��node����}����TP
		              while(node4 != NULL){
			          n4 = node4->GetValue()->GetNo();                                                                                                          //j'
			  
			                      char var_name[60];                                                                                                            //�ŧi�ܼƦW�٪��O���� ##
			                      sprintf(var_name,"V_%d_%d_%d",n2+1,n+1,n4+1);                                                                                 //��V_%d_%d_%d�g�Jvar_name�A�æL�Xj,i,j'
		 	                      bilvars.add(IloNumVar(env, BilVar_LB[ n2*M*N + n*N + n4 ], BilVar_UB[ n2*M*N + n*N + n4 ], var_name));                        //�NV�[�Jcplex���Ҥ�
			
			   
			     node4  = node4  ->GetNext();	
		        }
			 
		    node = node -> GetNext();		
		   }
		
		node2 = node2 -> GetNext();	
	   }
			 


//Add K=Poweri'j'*SINRj*U

	node2  = this-> GetNode_T();                                                                                                                                        //��node����}����TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                                                                                                                                //j
			//cout<< "n2=" << node2->GetValue()->GetNo() <<endl;
 		   node3 = this-> GetNode_B();                                                                                                                                  //��node����}����BS
			   while(node3 != NULL){
			   n3 = node3->GetValue()->GetNo();                                                                                                                         //i'
				//cout<< "n3=" << node3->GetValue()->GetNo() <<endl;
                  node4 = this-> GetNode_T();                                                                                                                           //��node����}����TP
		              while(node4 != NULL){
			          n4 = node4->GetValue()->GetNo();                                                                                                                  //j'

                                char var_name[60];                                                                                                                      //�ŧi�ܼƦW�٪��O���� ##
 			                    sprintf(var_name,"K_%d_%d_%d",n2+1,n3+1,n4+1 );                                                                                         //��K_%d_%d_%d�g�Jvar_name�A�æL�Xj,i',j'
	 		                    bilvars.add(IloNumVar(env,BilVar_LB[ N*M*N + n2*M*N + n3*N + n4 ], BilVar_UB[ N*M*N + n2*M*N + n3*N + n4 ], var_name));                 //�NK�[�Jcplex���Ҥ�
			                 
								

			  node4  = node4  ->GetNext();	
		     }     
			 
		   node3 = node3 -> GetNext();		
		  }
		
		node2 = node2 -> GetNext();	
	   } 
			 


 //Add ��V
			 
	node2  = this-> GetNode_T();                                 //��node����}����TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                         //j
	 		 
	node = this-> GetNode_B();                                   //��node����}����BS
		while(node != NULL){
		n = node->GetValue()->GetNo();                           //i
		 
	node4 = this-> GetNode_T();                                  //��node����}����TP
		while(node4 != NULL){
		n4 = node4->GetValue()->GetNo();                         //j'
		 
		for(int j=0 ; j<par ; j++)
		{

 				char var_name[60];
				sprintf(var_name,"��V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);                                                                                                         //��V�W�٤@�˪`�N ##
	 			bilvars.add(IloNumVar(env, BilVar_LB[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ], BilVar_UB[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ], var_name));     //�N��V�[�Jcplex���Ҥ�
			    
				//cout << "Sijn: BilVar_LB["<<i<<"*"<<M<<"*"<<N<<"*"<<par<<"+"<< i<<"*"<<N<<"*"<<par<<"+"<<n2<<"*"<<par<<"+"<<j<<"]="<< 4*M*N*M*N +n*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j<<endl;
	    }
	
			 node4 = node4 -> GetNext();	
		    }
			 
 
		     node = node -> GetNext();		
		    }
			 
 
			 node2 = node2 -> GetNext();	
		    }
			 
		 
			
//Add ��K
			 
	node2  = this-> GetNode_T();                                  //��node����}����TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                          //j
	 		 
	node3 = this-> GetNode_B();                                   //��node����}����BS
		while(node3 != NULL){
		n3 = node3->GetValue()->GetNo();                          //i'
		 

	node4 = this-> GetNode_T();                                   //��node����}����TP
		while(node4 != NULL){
		n4 = node4->GetValue()->GetNo();                          //j'
		 
		for(int j=0 ; j<par ; j++)
		{

 				char var_name[60];
				sprintf(var_name,"��K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                                                                                                                                 //��K�W�٤@�˪`�N ##
	 			bilvars.add(IloNumVar(env, BilVar_LB[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ], BilVar_UB[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ], var_name));    //�N��K�[�Jcplex���Ҥ�
			    
				//cout << "Sijn: BilVar_LB ="<<4*M*N*M*N + N*M*M*N*par + n*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j  <<endl;

		}

			 node4 = node4 -> GetNext();	
		    }
			 
		   
		     node3 = node3 -> GetNext();		
		    }
			 
 
			 node2 = node2 -> GetNext();	
		    }	 
			


			 
 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^---�ܼƫŧi����--------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 	//Object function ##############

 
	for(int i= 2*M*N + 2*N + M; i< 2*M*N + 3*N + M ; i++){
			Vars.add(vars[i]);             //Zj*u���O�����m
			Coes.add(Rate*k );          //R*b
		 	} 

 
		model.add(IloMaximize(env, IloScalProd(Coes,Vars)));

		Coes.clear();
		Vars.clear(); 

	  
//----------------------------------�ؼШ禡�ŧi����-------------------------------	

//----------------------------------------[Ij][Ji]�x�}�ŧi------------------------------------------

	double d;
	double Ij[200][200];																//��j��TP�i�H�s���Ҧ�BS
	double Ijcounter[200]
		;
	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			d = sqrt(pow(node2->GetValue()->GetX1() - node->GetValue()->GetX2(), 2) + pow(node2->GetValue()->GetY1() - node->GetValue()->GetY2(), 2) + pow((node2->GetValue()->GetZ1()) - node->GetValue()->GetZ2(), 2));


			if (d <= dmax)
			{

				Ij[n2][n] = n + 1;														//�Y�Z���bdmax�� �h��JBS�s��

				Ijcounter[n2] = Ijcounter[n2] + 1;

				cout << "Ij" << n2 << n << "=" << Ij[n2][n] << endl;

				cout << "Ijcounter" << n2 << "=" << Ijcounter[n2] << endl;


			}


			node2 = node2->GetNext();

		}
		node = node->GetNext();
	}


	//----------------------------------------------------------------------------------

	double Ji[200][200];
	double Jicounter[200];																			//��i��BS�i�H�s���Ҧ�TP

	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();
		//cout<< "n="<<node->GetValue()->GetNo()<<endl;

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();
			//cout<< "n2="<<node2->GetValue()->GetNo()<<endl;

			d = sqrt(pow(node2->GetValue()->GetX1() - node->GetValue()->GetX2(), 2) + pow(node2->GetValue()->GetY1() - node->GetValue()->GetY2(), 2) + pow((node2->GetValue()->GetZ1()) - node->GetValue()->GetZ2(), 2));
			//TP�PBS���Z��

			if (d <= dmax)
			{

				Ji[n][n2] = n2 + 1;																	//�Y�Z���bdmax�� �h��JTP�s��

				Jicounter[n] = Jicounter[n] + 1;

				cout << "Ji" << n << n2 << "=" << Ji[n][n2] << endl;

				cout << "Jicounter" << n << "=" << Jicounter[n] << endl;


			}


			node2 = node2->GetNext();

		}
		node = node->GetNext();
	}

	//----------------------------------------------------------------------------------
	//�x�sgain�x�}

	double gain[1000][1000];
	double Hnorm, d2, Gain2, Hnorm2;

	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();
		//cout<< "n="<<node->GetValue()->GetNo()<<endl;

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();

			d = sqrt(pow(node2->GetValue()->GetX1() - node->GetValue()->GetX2(), 2) + pow(node2->GetValue()->GetY1() - node->GetValue()->GetY2(), 2) + pow((node2->GetValue()->GetZ1()) - node->GetValue()->GetZ2(), 2));




			Hnorm = (137.4 + 35.2*(log10((d / 1000.0)))) / 10.0;//Hnorm=Lij/10   ;���|�l��Lij   ;d��쬰����
			gain[n][n2] = 1 / (pow(10.0, Hnorm));//�q�D�W�q

												 /*  Hnorm = (40*(log10(d/1000.0))+30*(log10(fa))+49)/10.0;




												 gain[n][n2] = 1/ (pow(10,Hnorm)) ; */

			cout << "gain" << n << n2 << "=" << gain[n][n2] << endl;
			cout << "d" << n << n2 << "=" << d << endl;

			node2 = node2->GetNext();

		}
		node = node->GetNext();
	}

	//---------------------------------------���-------------------------------------------

	//[A]���ƳW���O�s�ܼ�u��h�X�����
	node = this->GetNode_B();                 //��node����}����BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();          //i

		node2 = this->GetNode_T();                //��node����}����TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();        //j

			if (n == ((Ij[n2][n]) - 1))
			{
				if (n2 == ((Ji[n][n2]) - 1))
				{
					Coes.add(beta);
					Vars.add(vars[n*N + n2]);   //�s�WPower*u���ܼ�  
				}                                   //Coes�PVars�}�C���Ȧs�Y�ƻP�ܼ� �ϥΧ��M�� 
			}
			node2 = node2->GetNext();
		}                                           //beta*�U�UPij*u

		if (n == ((Ij[n2][n]) - 1))
		{
			Coes.add(Ps);                                        //Ps
			Vars.add(vars[2 * M*N + 2 * N + n]);                  //y*u
		}


		node = node->GetNext();
	}                                           //�UPs*yi*u

	char cst_name[60];
	sprintf(cst_name, "u");//��u�g�Jcst_name
	IloConstraint cst;//�Ыح��cst
	cst = IloScalProd(Coes, Vars) == 1;//���A
	cst.setName(cst_name);//�]�߭���W�� u
	model.add(cst); //�N����[�J�ҫ���
	Coes.clear();//�M��Coes�}�C
	Vars.clear();//�M��Vars�}�C


				 //[1]BS�S���}�ҡA�hTP�L�k�P��BS�i��s��
	node = this->GetNode_B();                          //��node����}����BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();                    //i

		node2 = this->GetNode_T();                      //��node����}����TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();              //j


			Coes.add(1);
			Vars.add(vars[M*N + 2 * N + n * N + n2]);   //�s�Wx*u���ܼ�
			Coes.add(-1);
			Vars.add(vars[2 * M*N + 2 * N + n]);        //�s�Wy*u���ܼ�



			sprintf(cst_name, "BS_%d_%d", n + 1, n2 + 1); //��BS_%d_%d�g�Jcst_name
			IloConstraint cst;                    //�Ыح��cst
			cst = IloScalProd(Coes, Vars) <= 0;     //���1 xij*u - yi*u <= 0
			cst.setName(cst_name);                //�]�߭���W��BS_%d_&d
			model.add(cst);                       //�N����[�J�ҫ���
			Coes.clear();                         //�M��Coes�}�C
			Vars.clear();                         //�M��Vars�}�C

			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}
	//�@M*N�����


	//[2]�ڭ̰��]��a�x���̤j�\�v�W����Pmax�P�ɡA��TP�S���s����Y��BS�ɡA��BS�]���|���t�U��ǿ�\�v����TP�C

	node = this->GetNode_B();                           //��node����}����BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();                    //i

		node2 = this->GetNode_T();                      //��node����}����TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();              //j

			Coes.add(-1.0);
			Vars.add(vars[n*N + n2]);           //�s�WPower*u���ܼ�

			Coes.add(Var_UB[n*N + n2]);
			Vars.add(vars[M*N + 2 * N + n * N + n2]);  //�s�Wx*u���ܼ�


			char cst_name[60];
			sprintf(cst_name, "Power_%d.%d", n + 1, n2 + 1); //��Power_%d.%d�g�Jcst_name
			IloConstraint cst;                    //�Ыح��cst
			cst = IloScalProd(Coes, Vars) >= 0;    //�ͦ��@�����
			cst.setName(cst_name);                //�]�߭���W��Power_%d
			model.add(cst);                       //�N����[�J�ҫ���
			Coes.clear();                         //�M��Coes�}�C
			Vars.clear();                         //�M��Vars�}�C

			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}
	//�@M*N�����


	//[3]BS�]�����O�b�}�Ҫ����p�U�A�~�i�H�A�ȩҳs����TP�C�BBS���t���Ҧ�TP���\�v�`�X�]����W�LPmax

	node = this->GetNode_B();                           //��node����}����BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();                     //i

		node2 = this->GetNode_T();                       //��node����}����TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();               //j

			Coes.add(1);
			Vars.add(vars[M*N + 2 * N + n * N + n2]);      //�s�WPower*u���ܼ�

			node2 = node2->GetNext();
		}
		Coes.add(-bsmaxpower);                           //BS��Max Power
		Vars.add(vars[2 * M*N + 2 * N + n]);           //�s�Wy*u���ܼ�  


		char cst_name[60];
		sprintf(cst_name, "SumPower_%d", n + 1); //��SumPower_%d�g�Jcst_name
		IloConstraint cst;                     //�Ыح��cst
		cst = IloScalProd(Coes, Vars) <= 0;      //�ͦ��@�����
		cst.setName(cst_name);                 //�]�߭���W��SumPower_%d
		model.add(cst);                        //�N����[�J�ҫ���
		Coes.clear();                          //�M��Coes�}�C
		Vars.clear();                          //�M��Vars�}�C

		node = node->GetNext();
	}
	//�@M�����


	//[4.1]�e
	//�YEj=0(���i�J�A��) �A�Ȩ��a�x�ƥج�0
	//�YEj=1�h�̤֥i�s1�Ӱ�a�x �̦h�i�H�sIj�x��a�x

	node2 = this->GetNode_T();                          //��node����}����TP
	while (node2 != NULL) {
		n2 = node2->GetValue()->GetNo();              //j

		node = this->GetNode_B();                       //��node����}����BS
		while (node != NULL) {
			n = node->GetValue()->GetNo();                //i
			if ((n == ((Ij[n2][n]) - 1)))
			{
				Coes.add(1.0);
				Vars.add(vars[M*N + 2 * N + n * N + n2]);   //�s�Wx*u���ܼ�
			}

			node = node->GetNext(); //�Uxij*u
		}

		Coes.add(-1.0);
		Vars.add(vars[2 * M*N + 2 * N + M + n2]);     //�s�WZj*u���ܼ�


		char cst_name[60];
		sprintf(cst_name, "TP-1-0_%d", n2 + 1);   //��TP_%d�g�Jcst_name
		IloConstraint cst;                    //�Ыح��cst
		cst = IloScalProd(Coes, Vars) >= 0;     //�ͦ��@����� �Uxij*u - zj*u >=0
		cst.setName(cst_name);                //�]�߭���W��TP_%d
		model.add(cst);                       //�N����[�J�ҫ���
		Coes.clear();                         //�M��Coes�}�C
		Vars.clear();	                      //�M��Vars�}�C

		node2 = node2->GetNext();
	}


	//[4.2]��

	node2 = this->GetNode_T();                          //��node����}����TP
	while (node2 != NULL) {
		n2 = node2->GetValue()->GetNo();               //j

		node = this->GetNode_B();                        //��node����}����BS
		while (node != NULL) {
			n = node->GetValue()->GetNo();                 //i
			if ((n == ((Ij[n2][n]) - 1)))
			{
				Coes.add(1.0);
				Vars.add(vars[M*N + 2 * N + n * N + n2]);    //�s�Wx*u���ܼ�
			}

			node = node->GetNext(); //�Uxij*u
		}

		Coes.add(-(Ijcounter[n2]));
		Vars.add(vars[2 * M*N + 2 * N + M + n2]);      //�s�WZj*u���ܼ�


		char cst_name[60];
		sprintf(cst_name, "TP-1-1_%d", n2 + 1);    //��TP_%d�g�Jcst_name
		IloConstraint cst;                     //�Ыح��cst
		cst = IloScalProd(Coes, Vars) <= 0;      //�ͦ��@����� �Uxij*u - Ij*Zj*u <= 0
		cst.setName(cst_name);                 //�]�߭���W��TP_%d
		model.add(cst);                        //�N����[�J�ҫ���
		Coes.clear();                          //�M��Coes�}�C
		Vars.clear();	                       //�M��Vars�}�C

		node2 = node2->GetNext();
	}											//�@N�����
  


//[5]���F�����Τ��t�Ϊ��ϥκ��N�סA�b��Ӧ�ʺ����л\�d�򤺡A�ܤ֭n�A�� mu%�H�W��TP�A�ȼƶq�A�䤤�Amu���w���`�ơA�N��̧CService Rate ���n�D�C 

node2  = this-> GetNode_T();                        //��node����}����TP
	  while (node2 != NULL){  
	  n2 = node2->GetValue()->GetNo();              //j
	       
				Coes.add(1);
				Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);    //�s�WZj*u���ܼ�					

			node2  = node2 ->GetNext();	
		   }	
		
			Coes.add(-(N*mu));
			Vars.add(vars[ 2*M*N + 3*N + M ]);      //�s�Wu���ܼ�

 
			char sct_name[60];
			sprintf(sct_name, "outtage ");          //��outtage �g�Jsct_name
			IloConstraint sct;
			sct = IloScalProd(Coes,Vars) >= 0 ;     //�ͦ��@�����
			sct.setName(sct_name);                  //�]�߭���W��outtage 
		    model.add(sct);                         //�N����[�J�ҫ���
			Coes.clear();                           //�M��Coes�}�C
			Vars.clear();                           //�M��Vars�}�C


//[6]Distance �䤤dmax��BS(i)�MTP(j)�����̤j�i�s�u�Z���A�O�Ӥw���`�ơC 

			double Gain;

	node  = this-> GetNode_B();                            //��node����}����BS
		while(node != NULL){
		n = node->GetValue()->GetNo();                      //i

		node2  = this-> GetNode_T();                        //��node����}����TP
			while (node2 != NULL){  
				n2 = node2->GetValue()->GetNo();            //j
				
				d = sqrt(pow(node2->GetValue()->GetX1()-node ->GetValue()->GetX2(),2)+pow(node2->GetValue()->GetY1()-node ->GetValue()->GetY2(),2)+pow((node2->GetValue()->GetZ1())-node ->GetValue()->GetZ2(),2));				//�Z��������}�ڸ�	
							
				Coes.add(d);
				Vars.add(vars[ M*N + 2*N + n*N + n2 ]);      //x*u

				Coes.add(-dmax);
				Vars.add(vars[ 2*M*N + 3*N + M ]);           //�s�Wu���ܼ�

					char var_name[60];
					sprintf(cst_name, "Distance_%d_%d", n+1, n2+1);     //��Distance_%d.%d�g�Jcst_name
					IloConstraint cst;                                  //�Ыح��cst
					cst = IloScalProd(Coes,Vars) <=0;                   //�ͦ��@�����
					cst.setName(cst_name);                              //�]�߭���W��Distance_%d.%d
					model.add(cst);                                     //�N����[�J�ҫ���
					Coes.clear();                                       //�M��Coes�}�C
					Vars.clear();                                       //�M��Vars�}�C
			
				node2  = node2  ->GetNext();
                 	
				}
			node  = node  ->GetNext();
			} 
		 
 			//�`�@M*N�����

 //[7]   //d2 �|���|�]�n�j  //���o��
//double d2, Gain2, Hnorm2,nfloor2;


node2  = this-> GetNode_T();
		while (node2 != NULL){  
	n2 = node2->GetValue()->GetNo(); 
 		//cout << "n2 =" <<  node2->GetValue()->GetNo()<<endl;

		node  = this-> GetNode_B();
			while(node != NULL){
			n = node->GetValue()->GetNo();
			//cout << "n =" <<  node->GetValue()->GetNo()<<endl;

			node4  = this-> GetNode_T();
				while (node4 != NULL){  
				n4 = node4->GetValue()->GetNo();  
               					//cout << "n4 =" <<  node4->GetValue()->GetNo()<<endl;

				node3  = this-> GetNode_B();
					while (node3 != NULL){
					n3 = node3->GetValue()->GetNo();
							//cout << "n3 =" <<  node3->GetValue()->GetNo()<<endl;

								
					d = sqrt(pow(node2->GetValue()->GetX1()-node ->GetValue()->GetX2(),2)+pow(node2->GetValue()->GetY1()-node ->GetValue()->GetY2(),2)+pow((node2->GetValue()->GetZ1())-node ->GetValue()->GetZ2(),2));							//�Z��������}�ڸ�



						/*	Hnorm= pow (10.0, 3.25*(log10(d/1000.0)) );
							Gain = 1.0/ (Keq*Hnorm) ;*/

	
					/*    Hnorm = (40*(log10(d/1000.0))+30*(log10(2300.0))+49)/10.0;
							  Gain = 1/ (pow(10,Hnorm)) ; */


							Hnorm = (137.4 + 35.2*(log10((d/1000.0))))/10.0;				//Lij(���|�l��)/10	
							Gain = 1.0/ (pow(10.0,Hnorm)) ;									//clannel gain

						/*	Hnorm = (143.503 + 38.35*(log10(d/1000.0)))/10.0;
							Gain = 1/ (pow(10,Hnorm)) ; */

		



						   d2 = sqrt(pow(node2->GetValue()->GetX1()-node3 ->GetValue()->GetX2(),2)+pow(node2->GetValue()->GetY1()-node3 ->GetValue()->GetY2(),2)+pow((node2->GetValue()->GetZ1())-node3 ->GetValue()->GetZ2(),2));				//�Z��������}�ڸ�



						/*		Hnorm2= pow (10.0, 3.25*(log10(d2/1000.0)) );
								Gain2 = 1.0/(Keq*Hnorm2) ;*/

			
						/*		Hnorm2 = (40*(log10(d2/1000.0))+30*(log10(2300.0))+49.0)/10.0;
			

			

								Gain2 = 1/ (pow(10,Hnorm2)) ;*/
								Hnorm2 = (137.4 + 35.2*(log10((d2/1000.0))))/10.0;			//�P�W
								Gain2 = 1.0/ (pow(10.0,Hnorm2)) ;

							/*  Hnorm2 = (143.503 + 38.35*(log10(d2/1000.0)))/10.0;
								Gain2 = 1/ (pow(10,Hnorm2)) ; */



								 /*        if(Ijcounter[n2]==3){
										   if (n!=n3 && n2!=n4 ){ //K
										   if((n!=((Ij[n2][n])-1))){
										   if ((n3!=((Ij[n2][n3])-1))&&(n4==(((Ji[n3][n4])-1)))){
		        
											Coes.add( Gain2/((2*pow(10.0,(-14)))) );
											Vars.add (bilvars[  M*N*N  + n2*M*N + n3*N + n4 ]);

											cout << "Gain2 =" <<  Gain2 <<endl;
											cout << "bilvars[  M*N*N  + n2*M*N + n3*N + n4 ] =" << bilvars[  M*N*N  + n2*M*N + n3*N + n4 ] <<endl;
											//Vars.add (bilvars[  M*N*M*N + n*N*M*N + n2*M*N + n3*N + n4 ]);
											//cout << "Gain2/Gain =" <<  Gain2/Gain <<endl;
											}}}}*/

									// 	    if(Ijcounter[n2]>=2){
											if (n==n3 && n2!=n4 ){												//i=i' j!=j'
									//		if((n!=((Ij[n2][n])-1))){
											 if ((n!=((Ij[n2][n])-1))&&(n4==(((Ji[n3][n4])-1)))){				//Ij�O��j��TP�i�H�s�W���Ҧ�BS , Ji�O��i��BS�i�H�s�W���Ҧ�TP	
																												//J���bI���л\�d�򤺡AJ'�bI���л\�d��
		        
												Coes.add( Gain2/((4.11*pow(10.0,(-14.0)))) );
												Vars.add (bilvars[  N*M*N + n2*M*N + n3*N + n4 ]);				//K

												cout << "Gain2 =" <<  Gain2 <<endl;
												cout << "bilvars[  N*M*N + n2*M*N + n3*N + n4 ] =" << bilvars[  N*M*N + n2*M*N + n3*N + n4 ] <<endl;
												//Vars.add (bilvars[  M*N*M*N + n*N*M*N + n2*M*N + n3*N + n4 ]);
												//cout << "Gain2/Gain =" <<  Gain2/Gain <<endl;
										   }}

						   
									/*     if(Ijcounter[n2]==1){
										   if (n!=n3 && n2!=n4 ){ //K
										   if((n==((Ij[n2][n])-1))){
										   if ((n3!=((Ij[n2][n3])-1))&&(n4==(((Ji[n3][n4])-1)))){
		        
												Coes.add( Gain2/((2*pow(10.0,(-14)))) );
												Vars.add (bilvars[  M*N*N  + n2*M*N + n3*N + n4 ]);

												cout << "Gain2 =" <<  Gain2 <<endl;
												cout << "bilvars[  M*N*N  + n2*M*N + n3*N + n4 ] =" << bilvars[  M*N*N  + n2*M*N + n3*N + n4 ] <<endl;
												//Vars.add (bilvars[  M*N*M*N + n*N*M*N + n2*M*N + n3*N + n4 ]);
												//cout << "Gain2/Gain =" <<  Gain2/Gain <<endl;
										   }}}}*/
      			
  		
										 if (n2!=n4 && n==n3 ){    												//i=i' j!=j'
										 if((n==((Ij[n2][n])-1))&&(n4==((Ji[n][n4])-1))){						//J�bI���л\�d�򤺡AJ'�bI���л\�d��
									//	 if ((n==(Ij[n2][n]-1))&& (n==n3) ){  //Z  	
									//	 if((n2!=n4)&&(n4==(Ji[n][n4]+1))){
									//   Coes.add( Alpha);
									//	 Coes.add(gain[n][n2]);
										 Coes.add( (Gain/(4.11*pow(10.0,(-14.0))) ));							//(4.11*pow(10.0,(-14.0) = Pn
										 Vars.add (bilvars[ n2*M*N + n*N + n4 ]);								//V
										 cout << "Alpha*Gain =" <<  Alpha*Gain <<endl;
										 cout << "bilvars[ n2*M*N + n*N + n4 ] =" << bilvars[ n2*M*N + n*N + n4 ] <<endl;
										// Vars.add (bilvars[ n*N*M*N + n2*M*N + n3*N + n4 ]); 
											 }
										 }
	 

											node3 = node3 -> GetNext();
					}     

										node4 = node4-> GetNext();
				}   


				if ((n==((Ij[n2][n])-1))){																		//�bi���л\�d��
						

					Coes.add(-(Gain/(4.11*pow(10.0,(-14.0))) ));												//(4.11*pow(10.0,(-14.0) = Pn
					Vars.add(vars[ n*N + n2 ]);																	//Power*U
					cout << "Gain =" << Gain <<endl;
					cout << "vars[ n*N + n2 ] =" <<  vars[ n*N + n2 ] <<endl;
				}
 
			 node = node -> GetNext();
			}
	   
    	
    // if ((n==((Ij[n2][n])-1))&&(n2==((Ji[n][n2])-1))){		
	//   Coes.add( Pn/Gain   );
	//Coes.add( 1 );
	//   Coes.add( Pn/gain[0][0] );
	//   Coes.add( Pn/(4.04182*pow(10.0,(-14.0))) );1.0*pow(10.0,(-10.0))
	    Coes.add(Pn/(4.11*pow(10.0,(-14.0))) );																	//Pn = 4.11*pow(10.0,(-14.0)	//�w�q�b�W��define
		Vars.add(vars[  M*N + n2 ]);																			//SINR*u

	//   cout << "(Pn/(2*pow(10.0,(-14))) )=" << (Pn/(2*pow(10.0,(-14))) ) <<endl;
	//Vars.add(vars[  M*N + n*N + n2 ]);   //sinr
		cout << "SINR_D_%d" << n2+1 <<endl;
		sprintf(cst_name, "SINR_D_%d",n2+1);
		IloConstraint cst;
		cst = IloScalProd(Coes,Vars) == 0;																		//���
		cst.setName(cst_name);
		model.add(cst);
		Vars.clear();
		Coes.clear();

		



	node2 = node2-> GetNext();
}  

//[8] Flow �ھ�Shannon capacity �A�ڭ̥i�H�p��BS i�MTP j ���q�D�e�q�C

   node2  = this-> GetNode_T();												//��node����}����TP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo();								//j
 
				Coes.add(-W / (log(2.0)));									//w/ln(2)
				Vars.add(vars[ M*N + N + n2 ]);								//H*u���O�����m

				Coes.add(Rate);												//R
				Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);						//Zj*u

				char cst_name[60];
				sprintf(cst_name, "Capacity_%d ", n2+1  );             
				IloConstraint cst;                                      //�Ыح��cst
				cst = IloScalProd(Coes,Vars) <= 0;                      //�ͦ��@�����
				cst.setName(cst_name);                                  
				model.add(cst);                                         //�N����[�J�ҫ���
			    Coes.clear();                                           //�M��Coes�}�C
                Vars.clear();                                           //�M��Vars�}�C
				  
	    	node2 = node2-> GetNext();
	      }   
	   


//-------------------------------------------------------Biliner Tearm : V--------------------------------------------

//[9.1] �fv*u
 
node  = this-> GetNode_B();                                             //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                      //i
	 
    node4  = this-> GetNode_T();                                        //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                //j'
 		 
				for(int j=0 ; j<par ; j++)
				{
					Coes.add(1);
				    Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //�s�W�fv*u���ܼ�
                } 

				Coes.add(-1);
				Vars.add(vars[ 2*M*N + 3*N + M ]);                                  //�s�Wu���ܼ�
				
				char cst_name[60];
				sprintf(cst_name, "Sum�fv_%d_%d ",n+1, n4+1 );          //��Sum�fv_%d.%d�g�Jcst_name
				IloConstraint cst;                                      //�Ыح��cst
				cst = IloScalProd(Coes,Vars) ==0;                       //�ͦ��@�����
				cst.setName(cst_name);                                  //�]�߭���W��
				model.add(cst);                                         //�N����[�J�ҫ���
			    Coes.clear();                                           //�M��Coes�}�C
                Vars.clear();                                           //�M��Vars�}�C
				  
		  
	    	node4 = node4-> GetNext();
	       }   
	   
    	    node = node -> GetNext();                                       
	       }
           //�`�@M*N�����
 
//[9.2.1]
			 
node2  = this-> GetNode_T();                                                       //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                               //j
	 
 	   node = this-> GetNode_B();                                                  //��node����}����BS
		   while(node != NULL){
		   n = node->GetValue()->GetNo();                                          //i
				
		      node4 = this-> GetNode_T();                                          //��node����}����TP
		            while(node4 != NULL){
			        n4 = node4->GetValue()->GetNo();                               //j'
			
			if (n2!=n4)                                                            // j!=j' 
			{                                                 
				for(int j=0 ; j<par ; j++)
			    {				 									
					Coes.add(-1);
	         		Vars.add (bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);					//�s�W��V*u par���ܼ� 
			    }  
		
				    Coes.add(1);
					Vars.add (vars[ M*N + n2 ]);                                                        //�s�WSINR*u���ܼ�
	
					sprintf(cst_name, "Sum_��V*u_%d_%d_%d",n2+1,n+1,n4+1);								//��Sum_��V*u_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) == 0;
					cst.setName(cst_name);
					model.add(cst);
				    Coes.clear();
     			    Vars.clear();	
			}				
		  
             node4 = node4 -> GetNext();	
		    }
			 
		     node = node -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			 
		   	//�`�@N*M*N�����
		    
 
//[9.2.2] SINR_��V
			 
node2  = this-> GetNode_T();                                                    //��node����}����TP
	 while(node2 != NULL){
	 n2 = node2->GetValue()->GetNo();                                           //j
	 
		node = this-> GetNode_B();                                              //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                                      //i
				 

		node4 = this-> GetNode_T();                                             //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();                                    //j'
		 
				    
			 if (n2!=n4)                                                        // j!=j' 
			 {                                               
				 for(int j=0 ; j<par ; j++)
				 {
 
					Coes.add(1);
	         		Vars.add (bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);					//�s�W��V*u par���ܼ�
		             
	          
					Coes.add(-(Var_UB [ M*N + n2 ] - Var_LB[ M*N + n2 ]));                              // (SINR_UB - SINR_LB)*u
					Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);						//�s�W�fv*u par���ܼ�
 
                    char cst_name[60];
			        sprintf(cst_name,  "��S_V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);							//�⡵S_V_%d_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0 ;
					cst.setName(cst_name);
					model.add(cst);
				    Coes.clear();
     			    Vars.clear();
 	            }		
			 } 
					  	              
			 node4  = node4  ->GetNext();	
		    }
			 
     	     node = node -> GetNext();		
		    }

	 
	         node2 = node2 -> GetNext();	
		    }
			 	             
            //�`�@N*M*N�����
 

//RLTV     
//[9.3.1] RLTV1		
 		 
 node2  = this-> GetNode_T();                                                 //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                          //j
	 	 
		node = this-> GetNode_B();                                            //��node����}����BS
		   while(node != NULL){
		   n = node->GetValue()->GetNo();                                     //i
				 
		       node4 = this-> GetNode_T();                                    //��node����}����TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();                           //j'

	 
		if (n2!=n4)     // j!=j' 
		{                                      
		 	for(int j=0 ; j<par ; j++)
			{
				    double a1 = ( ( j+1 ) / nsum );
				    double a2 =  ( j/nsum );
				    double a = ( pow( a1,ga ) - pow( a2,ga ) ) * (Var_UB[ n*N + n4 ] - Var_LB[ n*N + n4 ]);  


					Coes.add(-1*(Var_LB[ n*N + n2 ] + a*j ));                                                   //Power_LB + a*(np-1)
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);							//�s�W��V*u par���ܼ�
			}
	
					Coes.add(-1*Var_LB[ M*N + n2 ]);															//SINR_LB
					Vars.add(vars[ n*N + n4 ]);																	//�s�WPower*u���ܼ�

					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);													//�s�WV*u���ܼ�
                    
                    sprintf(cst_name, "V:RLT1_%d_%d_%d",n2+1,n+1,n4+1);											//��V:RLT1_%d_%d_%d�g�Jcst_name
                    IloConstraint cst;
				    cst = IloScalProd(Coes,Vars) >= 0;
				    cst.setName(cst_name);
				    model.add(cst);
				    Coes.clear();
     			    Vars.clear();	
	    		
		}		  
             node4 = node4 -> GetNext();	
		    }
			  
		     node = node -> GetNext();		
		    } 
		 
			 node2 = node2 -> GetNext();	
		    }
			 //�`�@N*M*N�����

  
//[9.3.2] RLTV2   

node2  = this-> GetNode_T();                                        //��node����}����TP
	 while(node2 != NULL){
     n2 = node2->GetValue()->GetNo();                               //j
	 		 
		node = this-> GetNode_B();                                  //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                          //i
				 
		       node4 = this-> GetNode_T();                          //��node����}����TP
		            while(node4 != NULL){
			        n4 = node4->GetValue()->GetNo();                //j'
	

		if (n2!=n4)        // j!=j' 
		{                                      
		 	for(int j=0 ; j<par ; j++)
			{
				    double a1 = (j+1)/nsum;
					double a2 = j/nsum;
					double a = ( pow( a1,ga ) - pow( a2,ga ) )  * (Var_UB[ n*N + n4 ] - Var_LB[ n*N + n4 ]); 


					Coes.add( -1 * ( Var_LB[ n*N + n4 ] + a*(j+1) ) );                                  //Power_LB + a*(np)
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);					//�s�W��V*u par���ܼ�
		    }
	
                  	Coes.add( -1*Var_LB[ M*N + n2 ]);                                                   //SINR_LB
					Vars.add(vars[ n*N + n4 ]);                                                         //�s�WPower*u���ܼ�
            
					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);                                            //�s�WV*u���ܼ�
				 
					char cst_name[60];
					sprintf(cst_name, "V:RLT2_%d_%d_%d", n2+1,n+1,n4+1);                                //��V:RLT2_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
     				Vars.clear();	
    		 		
	   } 
			 node4 = node4 -> GetNext();	
		    }
			 
 		   node = node -> GetNext();		
		  }
	
		node2 = node2 -> GetNext();	
	   } 
	   //�`�@N*M*N�����
 
//[9.3.3] RLTV3 
	 
node2  = this-> GetNode_T();                               //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                       //j
	 	 
		node = this-> GetNode_B();                         //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                 //i
				 
		       node4 = this-> GetNode_T();                 //��node����}����TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();        //j'
	 
		  if (n2!=n4)   // j!=j' 
		  {                  
			  for(int j=0 ; j<par ; j++)
			  {	
				    double a1 = ( ( j + 1 ) / nsum );
					double a2 = ( j / nsum );
					double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n*N + n4 ] - Var_LB[ n*N + n4 ]); 
					double Interval= Var_LB[ n*N + n2 ] + a*j;												//Interval = Power_LB + a*(np-1)


					Coes.add(-1*Interval);
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);						//�s�W��V*u par���ܼ�

					Coes.add ( Interval*(Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ]) );                        //Interval*(SINR_UB-SINR_LB)
					Vars.add (vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);                          //�s�W�fv*u par���ܼ� 
		 
         	  }	
		  					
					Coes.add(-1*Var_UB[ M*N + n2 ]);                                                        //SINR_UB
					Vars.add(vars[ n*N + n4 ]);                                                             //�s�WPower*u���ܼ�

					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);												//�s�WV*u���ܼ�
				    
					sprintf(cst_name, "V:RLT3_%d_%d_%d",n2+1,n+1,n4+1 );                                    //��V:RLT3_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
     				Vars.clear();	
		  }
			 
			 node4 = node4 -> GetNext();	
		    }
			 
 		   node = node -> GetNext();		
		  }
	
		node2 = node2 -> GetNext();	
	   }
	   //�`�@N*M*N�����
	
  
//[9.3.4] RLTV4
 		 
node2  = this-> GetNode_T();                               //��node����}����TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                       //j
	 
 	   node = this-> GetNode_B();                          //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                 //i
				 
		      node4 = this-> GetNode_T();                  //��node����}����TP
		          while(node4 != NULL){
			      n4 = node4->GetValue()->GetNo();         //j'


		  if (n2!=n4) // j!=j' 
		  {                 
			  for(int j=0 ; j<par ; j++)
			  {
				    double a1 = (j+1)/nsum;
					double a2 = j/nsum;
					double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n*N + n4 ] - Var_LB[ n*N + n4 ]);
					double Interval= Var_LB[ n*N + n2 ] + a*(j+1);												//Interval = Power_LB + a*(np)


					Coes.add(-1*Interval);
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);							//�s�W��V*u par���ܼ�

					Coes.add ( Interval*(Var_UB[ M*N + n2 ]- Var_LB[ M*N + n2 ]));								//Interval*(SINR_UB-SINR_LB)
					Vars.add ( vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);								//�s�W�fv*u par���ܼ�
		 
         	  }	
		  

					Coes.add(-1*Var_UB[ M*N + n2 ]);                                                            //SINR_UB
					Vars.add(vars[ n*N + n4 ]);                                                                 //�s�WPower(i,j')*u���ܼ�

					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);                                                    //�s�WV*u���ܼ�
				    
					sprintf(cst_name, "V:RLT4_%d_%d_%d",n2+1,n+1, n4+1 );										//��V:RLT4_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) >=0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
     				Vars.clear();	
		 }
			 
			 node4 = node4 -> GetNext();	
		    }
			 
 		   node = node -> GetNext();		
		  }
	
		node2 = node2 -> GetNext();	
	   }
	   //�`�@N*M*N����� 
	
    
//-------------------------------------------------------Bilinear Tearm :K---------------------------------------//

//[10.1] �fk*u
 
node3  = this-> GetNode_B();                                             //��node����}����BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();                                     //i'
	 
    node4  = this-> GetNode_T();                                         //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                 //j'
 		 
				for(int j=0 ; j<par ; j++)
				{
					Coes.add(1);
				    Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);   //�s�W�fk*u���ܼ�
                } 

				Coes.add(-1);
				Vars.add(vars[ 2*M*N + 3*N + M ]);                                             //�s�Wu���ܼ�
				
				char cst_name[60];
				sprintf(cst_name, "Sum�fk_%d_%d ",n3+1, n4+1 );         //��Sum�fk_%d.%d�g�Jcst_name
				IloConstraint cst;                                      //�Ыح��cst
				cst = IloScalProd(Coes,Vars) ==0.0;                     //�ͦ��@�����
				cst.setName(cst_name);                                  //�]�߭���W��Sum�f_%d_%d
				model.add(cst);                                         //�N����[�J�ҫ���
			    Coes.clear();                                           //�M��Coes�}�C
                Vars.clear();                                           //�M��Vars�}�C
				  
		  
	    	node4 = node4 -> GetNext();
	       }   
	   
    	    node3 = node3 -> GetNext();                                       
	       }
           //�`�@M*N�����


//[10.2.1]
			  
node2  = this-> GetNode_T();                                     //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                             //j
	 	 
		node3 = this-> GetNode_B();                              //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                     //i'
				 
		       node4 = this-> GetNode_T();                       //��node����}����TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();              //j'
			
			if (n2!=n4)        // j!=j' 
			{                       
				for(int j=0 ; j<par ; j++)
				{					
					
					Coes.add(-1);    
	         		Vars.add (bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);       //�s�W��K*u���ܼ�	  
		        
				}  
                
				    Coes.add(1);
					Vars.add (vars[ M*N + n2 ]);									//�s�WSINR*u���ܼ� 
                    	
					sprintf(cst_name, "Sum_��K*u_%d_%d_%d",n2+1,n3+1,n4+1);			//��Sum_��K*u_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) == 0;
					cst.setName(cst_name);
					model.add(cst);
				    Coes.clear();
     			    Vars.clear();	
			}			
		   
             node4 = node4 -> GetNext();	
		    }
			 
		     node3 = node3 -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N�����
		    

//[10.2.2] SINR_��K   
 		 
node2  = this-> GetNode_T();                             //��node����}����TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                     //j
	 	
		node3 = this-> GetNode_B();                      //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();             //i'
				 
		       node4 = this-> GetNode_T();               //��node����}����TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();      //j'
		 

		if (n2!=n4)      // j!=j'
		{                         
			for(int j=0 ; j<par ; j++)
			{
										
				    Coes.add(1);                           
	         		Vars.add (bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);			//�s�W��K*u���ܼ�
					                 
				    Coes.add(-(Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ]));
				    Vars.add( vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);				//�s�W�fk*u���ܼ�

					sprintf(cst_name,  "��S_K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                             //�⡵S_K_%d_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0 ;
					cst.setName(cst_name);
					model.add(cst);
				    Coes.clear();
     			    Vars.clear();
	
 	         }
			
	    } 
				
		 
             node4 = node4 -> GetNext();	
		    }
			 
		  
		     node3 = node3 -> GetNext();		
		    }
			 
			
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N�����
  
//RLTK
//[10.3.1] RLTK1  		
 	 
node2  = this-> GetNode_T();                                              //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                      //j
	 
 		 
	   node3 = this-> GetNode_B();                                        //��node����}����BS
		   while(node3 != NULL){
		   n3 = node3->GetValue()->GetNo();                               //i'
				 

		     node4 = this-> GetNode_T();                                  //��node����}����TP
		         while(node4 != NULL){
			     n4 = node4->GetValue()->GetNo();                         //j'
	 

		if (n2!=n4)    // j!=j'
		{                   

			for(int j=0 ; j<par ; j++)
			{
				    double a1 = (j+1)/nsum;
				    double a2 = j/nsum;
				    double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n3*N + n4 ] - Var_LB[ n3*N + n4 ]);		//��k1  *


					Coes.add( -1*(Var_LB[ n3*N + n2 ] + a*j) );                                         //Power_LB + a*(np-1)
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);		//�s�W��K*u���ܼ�

			}	
				    
			        Coes.add(-1*Var_LB[ M*N + n2 ]);													//SINR_LB
					Vars.add(vars[ n3*N + n4 ]);														//�s�WPower*u���ܼ�

					Coes.add(1);
  					Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);									//�s�WK*u���ܼ�
			 
                    sprintf(cst_name, "K:RLT1_%d_%d_%d", n2+1,n3+1,n4+1);                               //��K:RLT1_%d_%d_%d�g�Jcst_name
				    IloConstraint cst;
				    cst = IloScalProd(Coes,Vars) >= 0;
				    cst.setName(cst_name);
				    model.add(cst);
				    Coes.clear();
     			    Vars.clear();	
	    }		
				  
             node4 = node4 -> GetNext();	
		    }
			 
		     node3 = node3 -> GetNext();		
		    }
					
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N����� 
		   	  
   
//[10.3.2] RLTK2
			 
node2  = this-> GetNode_T();                                      //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                              //j
	 
		node3 = this-> GetNode_B();                               //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                      //i'
				 
               node4 = this-> GetNode_T();                        //��node����}����TP
	              while(node4 != NULL){
			      n4 = node4->GetValue()->GetNo();                //j'


		if (n2!=n4 )         // j!=j'
		{                                      
			for(int j=0 ; j<par ; j++)
			{
				    double a1 = ((j+1) / nsum);
					double a2 = (j / nsum);
					double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n3*N + n4 ] - Var_LB[ n3*N + n4 ]);				  

		            Coes.add( -1*(Var_LB[ n3*N + n2 ] + a*(j+1) ) );                                    //Power_LB + a*np
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j  ]);     //�s�W��K*u���ܼ�
		 
		    }	
	
			        Coes.add(-1*Var_LB[ M*N + n2]);														//SINR_LB
					Vars.add(vars[ n3*N + n4 ]);														//�s�WPower*u���ܼ�

					Coes.add(1);
  					Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);                                   //�s�WK*u���ܼ�		

                    sprintf(cst_name, "K:RLT2_%d_%d_%d",n2+1,n3+1,n4+1);								//��K:RLT2_%d_%d_%d�g�Jcst_name
			        IloConstraint cst;
				    cst = IloScalProd(Coes,Vars) <= 0 ;
				    cst.setName(cst_name);
				    model.add(cst);
				    Coes.clear();
     			    Vars.clear();	
						
		 }	  
             node4 = node4 -> GetNext();	
		    }
		
		     node3 = node3 -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N����� 
	

//[10.3.3] RLTK3 
		 
node2  = this-> GetNode_T();                                          //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                  //j
	 
 		 
		node3 = this-> GetNode_B();                                   //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                          //i'
				 

		       node4 = this-> GetNode_T();                            //��node����}����TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();                   //j'
	 

		  if (n2!=n4)       // j!=j'
		  {                        
			  for(int j=0 ; j<par ; j++)
			  {
				    double a1 = ((j+1) / nsum);
					double a2 = (j / nsum);
					double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n3*N + n4 ] - Var_LB[ n3*N + n4 ]); 

					double Interval= Var_LB[ n3*N + n2 ] + a*j;
			
					Coes.add(-1*Interval);
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);		//�s�W��K*u���ܼ�
			 

					Coes.add ( Interval*(Var_UB[ M*N + n2]- Var_LB[ M*N + n2]));						//Interval*(SINR_UB-SINR_LB)
					Vars.add ( vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);               //�s�W�fk*u���ܼ�
		 
              }	
						
				    Coes.add(-1*Var_UB[ M*N + n2 ]);													//SINR_UB
					Vars.add(vars[ n3*N + n4 ]);														//�s�WPower*u���ܼ�

					Coes.add(1);
	         		Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);									//�s�WK*u���ܼ�
				    
					sprintf(cst_name, "K:RLT3_%d_%d_%d", n2+1,n3+1,n4+1 );								//��K:RLT3_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
     				Vars.clear();	
		  }
			 
			 node4  = node4  ->GetNext();	
		    }
			 
 		     node3 = node3 -> GetNext();		
		    }
	
		     node2 = node2 -> GetNext();	
	        }
			//�`�@N*M*N����� 
	
  
//[10.3.4] RLTK4
			 
node2  = this-> GetNode_T();                                   //��node����}����TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                           //j
	 		 
		node3 = this-> GetNode_B();                            //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                   //i'
				 
		       node4 = this-> GetNode_T();                     //��node����}����TP
		          while(node4 != NULL){
			      n4 = node4->GetValue()->GetNo();             //j'
   

		  if (n2!=n4)      //j!=j'
		  {                      
			  for(int j=0 ; j<par ; j++)
			  {
				    double a1 = (j+1)/nsum;
					double a2 = j/nsum;
					double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n3*N + n4 ] - Var_LB[ n3*N + n4 ]);

					double Interval= Var_LB[ n3*N + n2 ] + a*(j+1);										//Interval = Power_LB + a*np

		
					Coes.add(-1*Interval);
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);		//�s�W��K*u���ܼ�
		 

					Coes.add ( Interval*(Var_UB[ M*N + n2 ]- Var_LB[ M*N + n2 ]));                      //Interval*(SINR_UB-SINR_LB)
					Vars.add ( vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);          //�s�W�fk*u���ܼ�
		 
         	  }	
		  				
				    Coes.add(-1*Var_UB[ M*N + n2]);														//SINR_UB
					Vars.add(vars[ n3*N + n4 ]);														//�s�WPower*u���ܼ�

					Coes.add(1);
	         		Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);									//�s�WK*u���ܼ�
				 
					sprintf(cst_name, "K:RLT4_%d_%d_%d",n2+1,n3+1, n4+1 );								//��K:RLT4_%d_%d_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) >= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
     				Vars.clear();	
		  }
			 
			node4  = node4  ->GetNext();	
		   }
			 
 			node3 = node3 -> GetNext();		
		   }
	
		    node2 = node2 -> GetNext();	
	       }
           //�`�@N*M*N�����
 

//[12.2 - 12.5] Capacity (�]�t������ƶ��ءA�n�O�o�ݨD�Ѥ覡)
 
    float HUB, HLB, slmMinus,R , D, Q , Padd1;	
 

//[12.2]
			 
node2  = this-> GetNode_T();										//��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();								//j
	 
 		 
				
					 Padd1 = 1 + Var_LB[ M*N + n2 ];                //1+SINR_LB
					 R = log(Padd1)-1;
					 Q = R*Padd1 + 1;
					
					//cout<< "Pmax =" <<Pmax<<endl;

					Coes.add(Padd1);                                //1+SINR_LB
					Vars.add(vars[ M*N + N + n2 ]);					//�s�WH*u���ܼ�
					
					Coes.add(-1);
					Vars.add(vars[ M*N + n2 ]);						//�s�WSINR*u���ܼ�
					
					Coes.add(-Q);
					Vars.add(vars[ 2*M*N + 3*N + M ]);              //�s�Wu���ܼ�
										
                    
					char cst_name[60];
					sprintf(cst_name, "Capacity1_D_%d",n2+1);		//��Capacity1_D_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0 ;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();

			 node2 = node2 -> GetNext();	
		    }
			 
            
	        //�`�@N�����

//[12.3]
		 
node2  = this-> GetNode_T();             //��node����}����TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();     //j
	 
				    HUB = 1 + Var_UB[ M*N + n2 ];							//HUB=1+SINR_UB
				    HLB = 1 + Var_LB[ M*N + n2 ];							//HLB=1+SINR_LB
				    slmMinus = Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ];		//SINR_UB-SINR_LB
			
 
				    Padd1 = 1 + beta;										//beta���ؼЦ�				
					R = log(Padd1)-1;
					Q = R*Padd1 + 1;
 
					Coes.add(Padd1);
					Vars.add(vars[ M*N + N + n2 ]);							//�s�WH*u���ܼ�
					
					Coes.add(-1);
					Vars.add(vars[ M*N + n2 ]);								//�s�WSINR*u���ܼ�
					
					Coes.add(-Q);
					Vars.add(vars[ 2*M*N + 3*N + M ]);						//�s�Wu���ܼ�
				
                    char cst_name[60];
					sprintf(cst_name, "Capacity2_%d", n2+1);				//��Capacity2_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();
				
			 node2 = node2 -> GetNext();	
		    }
	
	        //�`�@N�����
					
//[12.4]
			 
node2  = this-> GetNode_T();												//��node����}����TP  
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();										//j

					Padd1 = 1 + Var_UB[ M*N + n2 ];							//Padd1=1+SINR_UB
					R = log(Padd1)-1;
					Q = R*Padd1 + 1;
										
					Coes.add(Padd1);
					Vars.add(vars[ M*N + N + n2 ]);							//�s�WH*u���ܼ�
					
					Coes.add(-1);
					Vars.add(vars[ M*N + n2 ]);								//�s�WSINR*u���ܼ�

					Coes.add(-Q);
					Vars.add(vars[ 2*M*N + 3*N + M ]);						//�s�Wu���ܼ�
 
                    char cst_name[60];
					sprintf(cst_name, "Capacity3_%d", n2+1);				//��Capacity3_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();
					
			 node2 = node2 -> GetNext();	
		    } 
			 
	        //�`�@N�����
					
//[12.5] 
			 
node2  = this-> GetNode_T();                 //��node����}����TP
	while(node2 != NULL){  
	n2 = node2->GetValue()->GetNo();         //j
			      
					slmMinus = Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ];         //SINR_UB-SINR_LB
		     		HUB = log(1 + Var_UB[ M*N + n2 ] );                         //ln(1+SINR_UB)
					HLB = log(1 + Var_LB[ M*N + n2 ] );                         //ln(1+SINR_LB)
					R =  Var_UB[ M*N + n2 ]*HLB;                                //SINR_UB*HLB
					D =  Var_LB[ M*N + n2 ]*HUB;                                //SINR_LB*HUB
            

					Coes.add(slmMinus);
					Vars.add(vars[ M*N + N + n2 ]);								//�s�WH*u���ܼ�

					Coes.add(HLB-HUB);
					Vars.add(vars[ M*N + n2 ]);									//�s�WSINR*u���ܼ�

					Coes.add(D-R);
					Vars.add(vars[ 2*M*N + 3*N + M ]);							//�s�Wu���ܼ�

                    char cst_name[60];
		            sprintf(cst_name, "Capacity4_%d", n2+1);					//��Capacity4_%d�g�Jcst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) >= 0 ;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();
				
		     node2 = node2 -> GetNext();	
		    }
			
            //�`�@N�����

//��u���
//[13.1] y*u

 node  = this-> GetNode_B();                   //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();             //i
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + n ]);       //�s�Wy*u���ܼ�

			Coes.add(-1);                      
			Vars.add(vars[ 2*M*N + 3*N + M ]);       //�s�Wu���ܼ�

			sprintf(cst_name, "y_u_%d",n+1);         //��y_u_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node  = node  -> GetNext();
	    }
	    //�`�@M����� 

//[13.2] x*u

 node  = this-> GetNode_B();               //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();         //i
     	 	    
		node2  = this-> GetNode_T();       //��node����}����TP
		while (node2 != NULL){  
		n2 = node2->GetValue()->GetNo();   //j 
		    	

			Coes.add(1);
			Vars.add(vars[ M*N + 2*N + n*N + n2 ]);        //�s�Wx*u���ܼ�

			Coes.add(-1); 
			Vars.add(vars[ 2*M*N + 3*N + M ]);             //�s�Wu���ܼ�

			sprintf(cst_name, "x_u_%d_%d",n+1,n2+1);       //��x_u_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <=0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		    node2 = node2 -> GetNext();
		   }
		    node  = node  -> GetNext();
	       }
	       //�`�@M*N�����

//[13.3] Zj*u

 node2  = this-> GetNode_T();                   //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();            //j
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);   //�s�WZj*u���ܼ�

			Coes.add(-1);                      
			Vars.add(vars[ 2*M*N + 3*N + M ]);        //�s�Wu���ܼ�
            
			sprintf(cst_name, "Zj_u_%d",n2+1);        //��Zj_u_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node2  = node2  -> GetNext();
	    }
	    //�`�@N�����
 
//[13.4.1] �fv*u

 node  = this-> GetNode_B();                      //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                //i
     	 	    
		node4  = this-> GetNode_T();              //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();          //j' 
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //�s�W�fv*u���ܼ�

			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);                              //�s�Wu���ܼ�
		
			sprintf(cst_name, "�fv_u_%d_%d_%d",n+1,n4+1,j+1);               //��fv_u_%d_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <=0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

	    }
					 
		   node4 = node4 -> GetNext();
		  }
		   node  = node  -> GetNext();
	      }
	      //�`�@M*N�����

//[13.4.2] �fk*u

 node3 = this-> GetNode_B();                      //��node����}����BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();              //i'
     	 	    
		node4  = this-> GetNode_T();              //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();          //j' 
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);   //�s�W�fk*u���ܼ�

			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);                                         //�s�Wu���ܼ�
		
			sprintf(cst_name, "�fk_u_%d_%d_%d",n3+1,n4+1,j+1);                         //��fk_u_%d_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <=0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

	    }
					 
		   node4 = node4 -> GetNext();
		  }
		   node3 = node3 -> GetNext();
	      }
	      //�`�@M*N�����


//[13.5] ��V*u
/*		 
node2  = this-> GetNode_T();                     //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();             //j
	 		 
	    node = this-> GetNode_B();               //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();       //i
				 
		node4 = this-> GetNode_T();              //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();     //j'
        
        
			for(int j=0 ; j<par ; j++)
			{
				Coes.add(1);
			    Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);   //�s�W��V*u���ܼ�

			    Coes.add(-1);
			    Vars.add(vars[ 2*M*N + 3*N + M ]);                                                //�s�Wu���ܼ�
		
			    sprintf(cst_name, "��V_u_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1 );                        //�⡵V_u_%d_%d_%d_%d�g�Jcst_name
 			    IloConstraint cst;
			    cst = IloScalProd(Coes,Vars) <=0;
			    cst.setName(cst_name);
			    model.add(cst); 
		   	    Coes.clear();
			    Vars.clear();	
			}	
		 
		     node4  = node4  ->GetNext();	
		    }
			 
		     node = node -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N�����
 
//[13.6] ��K*u
			 
node2  = this-> GetNode_T();                   //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();           //j
	 
		node3 = this-> GetNode_B();            //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();   //i'
				 
		node4 = this-> GetNode_T();            //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();   //j'
		 
        
			for(int j=0 ; j<par ; j++)
			{
				Coes.add(1);
			    Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);    //�s�W��K*u���ܼ�

			    Coes.add(-1);
		      	Vars.add(vars[ 2*M*N + 3*N + M ]);                                                  //�s�Wu���ܼ�
		
			    sprintf(cst_name, "��K_u_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1 );                         //�⡵K_u_%d_%d_%d_%d�g�Jcst_name
 			    IloConstraint cst;
			    cst = IloScalProd(Coes,Vars) <=0;
			    cst.setName(cst_name);
			    model.add(cst); 
		   	    Coes.clear();
			    Vars.clear();	
		    }
	    		 
		     node4 = node4 -> GetNext();	
		    }
			 
		     node3 = node3 -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N*par�����

*/
 
//[14] Big M-------------------------------------------------------------------------
 
//[14.1] y*u

 node  = this-> GetNode_B();                      //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                //i
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + n ]);          //�s�Wy*u���ܼ�

			Coes.add(-Big);                             //Big M
			Vars.add(binvars[ n ]);                     //�s�Wy���ܼ�
            
			sprintf(cst_name, "Big M_y_%d",n+1);        //��Big M_y_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

	    node  = node  -> GetNext();
	   }
	   //�`�@M�����
 
//[14.2] x*u

 node  = this-> GetNode_B();                       //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                 //i
     	 	    
		node2  = this-> GetNode_T();               //��node����}����TP
		while (node2 != NULL){  
		n2 = node2->GetValue()->GetNo();           //j
		    	

			Coes.add(1);
			Vars.add(vars[ M*N + 2*N + n*N + n2 ]);    //�s�Wx*u���ܼ�

			Coes.add(-Big);                            //Big M
			Vars.add(binvars[ M + n*N + n2 ]);         //�s�Wx���ܼ�

			sprintf(cst_name, "Big M_x_%d_%d",n+1,n2+1);   //��Big M_x_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <=0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	
					 
		    node2 = node2 -> GetNext();
		   }
		    node  = node  -> GetNext();
	       }
	       //�`�@M*N�����

//[14.3] Zj*u

 node2 = this-> GetNode_T();                      //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();              //j
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);       //�s�WZj*u���ܼ�

			Coes.add(-Big);                         //Big M
			Vars.add(binvars[ M + M*N + n2 ]);      //�s�WZj���ܼ�

			sprintf(cst_name, "Big M_Zj_%d",n2+1);  //��Big M_Zj_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

	    node2 = node2 -> GetNext();
	   }
	   //�`�@N�����
 
//[14.4.1] �fv*u

 node  = this-> GetNode_B();                        //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                  //i
     	 	    
		node4  = this-> GetNode_T();                //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();            //j'
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //�s�W�fv*u���ܼ�

			Coes.add(-Big);                                           //Big M
			Vars.add(binvars[ M + M*N + N + n*N*par + n4*par + j ]);  //�s�W�fv���ܼ�
		
			sprintf(cst_name, "Big M_�fv_%d_%d_%d",n+1,n4+1,j+1);     //��Big M_�fv_%d_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <=0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		}
					 
		   node4 = node4 -> GetNext();
		  }
		   node  = node  -> GetNext();
	      }
	      //�`�@M*N�����
  
//[14.4.2] �fk*u

 node3 = this-> GetNode_B();                        //��node����}����BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();                //i'
     	 	    
		node4  = this-> GetNode_T();                //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();            //j'
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);   //�s�W�fk*u���ܼ�

			Coes.add(-Big);                                                      //Big M
			Vars.add(binvars[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ]);  //�s�W�fk���ܼ�
		
			sprintf(cst_name, "Big M_�fk_%d_%d_%d",n3+1,n4+1,j+1);    //��Big M_�fk_%d_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <=0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		}
					 
		   node4 = node4 -> GetNext();
		  }
		   node3 = node3 -> GetNext();
	      }
	      //�`�@M*N�����
  

//[14.5] ��V*u  
/*		 
node2  = this-> GetNode_T();                       //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();               //j
 		 
		node = this-> GetNode_B();                 //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();         //i
				 
		node4 = this-> GetNode_T();                //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();       //j'
		   
		if(n2!=n4) // j!=j'
		{                            
			for(int j=0 ; j<par ; j++)
		    {
				Coes.add(1);
			    Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);   //�s�W��V*u���ܼ�

			    Coes.add(-Big);                                                                   //Big M
			    Vars.add(bilvars[ 2*N*M*N + 2*N*M*N*par + n2*M*N*par + n*N*par + n4*par + j ]);                 //�s�W

				sprintf(cst_name, "Big M_��V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);                     //��Big M_��V_%d_%d_%d_%d�g�Jcst_name
 			    IloConstraint cst; 
			    cst = IloScalProd(Coes,Vars) <=0;
			    cst.setName(cst_name);
			    model.add(cst); 
		   	    Coes.clear();
			    Vars.clear();	
			}

		}	

		     node4  = node4  ->GetNext();	
		    }
			 
		     node = node -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
	        //�`�@N*M*N�����
 

//[14.6] ��K*u 
	 
	node2  = this-> GetNode_T();              //��node����}����TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();      //j
 		 
		node3 = this-> GetNode_B();           //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();  //i'
				 
		node4 = this-> GetNode_T();           //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();  //j'

			if(n2!=n4)    // j!=j'
			{
				for(int j=0 ; j<par ; j++)
		        {
					Coes.add(1);
			        Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);      //�s�W��K*u���ܼ�

			        Coes.add(-Big);                                                                       //Big M
			        Vars.add(bilvars[ 2*N*M*N + 3*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);        //�s�W
 
					sprintf(cst_name, "Big M_��K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                        //��Big M_��K_%d_%d_%d_%d�g�Jcst_name
 			        IloConstraint cst;  
			        cst = IloScalProd(Coes,Vars) <=0;
			        cst.setName(cst_name);
			        model.add(cst); 
		   	        Coes.clear();
			        Vars.clear();	

		       }

			}	

		     node4 = node4 -> GetNext();	
		    }
			 
		     node3 = node3 -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N�����
		   

 */
//[15] Big M 
 
//[15.1] y*u

 node  = this-> GetNode_B();                       //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                 //i
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);           //�s�Wu���ܼ�

			Coes.add(Big);                         //Big M
			Vars.add(binvars[ n ]);                //�s�Wy���ܼ�
			
			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 2*N + n ]);     //�s�Wy*u���ܼ�
            
			sprintf(cst_name, "Mix M_y_%d",n+1);   //��Mix M_y_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node  = node  -> GetNext();
	   }
       //�`�@M�����


//[15.2] x*u

 node  = this-> GetNode_B();                     //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();               //i
     	 	    
		node2  = this-> GetNode_T();             //��node����}����TP
		while (node2 != NULL){  
		n2 = node2->GetValue()->GetNo();         //j
		    	

		    Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);          //�s�Wu���ܼ�

			Coes.add(Big);                        //Big M
			Vars.add(binvars[ M + n*N + n2 ]);    //�s�Wx���ܼ�
			
			Coes.add(-1);
     		Vars.add(vars[ M*N + 2*N + n*N + n2 ]);   //�s�Wx*u���ܼ�
 
			sprintf(cst_name, "Mix_M_x_%d_%d",n+1,n2+1);  //��Mix_M_x_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

					 
		    node2 = node2 -> GetNext();
		   }
		    node  = node  -> GetNext();
	       }
	       //�`�@M*N�����


//[15.3] Zj*u

 node2  = this-> GetNode_T();                        //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                 //j
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);             //�s�Wu���ܼ�

			Coes.add(Big);                           //Big M
			Vars.add(binvars[ M + M*N + n2 ]);       //�s�W��Zj�ܼ�
			
			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);        //�s�WZj*u���ܼ�
            
			sprintf(cst_name, "Mix M_Zj_%d",n2+1);   //��Mix M_Zj_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node2  = node2 -> GetNext();
	   }
       //�`�@N�����

 
//[15.4.1] �fv*u

 node  = this-> GetNode_B();                  //��node����}����BS
	while(node != NULL){
	n = node->GetValue()->GetNo();            //i
     	 	    
		node4  = this-> GetNode_T();          //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();      //j'
		    	
		for(int j=0 ; j<par ; j++)
		{

			Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);								//�s�Wu���ܼ�

			Coes.add(-1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //�s�W�fv*u���ܼ�
 
			Coes.add(Big);													//Big M
			Vars.add(binvars[ M + M*N + N + n*N*par + n4*par + j ]);		//�s�W�fv���ܼ�
		
			sprintf(cst_name, "Mix_M_�fv_%d_%d_%d",n+1,n4+1,j+1);			//��Mix_M_�fv_%d_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		}
					 
		    node4 = node4 -> GetNext();
		   }
		    node  = node  -> GetNext();
	       }
           //�`�@M*N�����


//[15.4.2] �fk*u

 node3  = this-> GetNode_B();                  //��node����}����BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();           //i'
     	 	    
		node4  = this-> GetNode_T();           //��node����}����TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();       //j'
		    	
		for(int j=0 ; j<par ; j++)
		{

			Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);      //�s�Wu���ܼ�

			Coes.add(-1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);    //�s�W�fk*u���ܼ�
 
			Coes.add(Big);                                                        //Big M
			Vars.add(binvars[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ]);   //�s�W�fk���ܼ�
		
			sprintf(cst_name, "Mix_M_�fk_%d_%d_%d",n3+1,n4+1,j+1);                //��Mix_M_�fk_%d_%d_%d�g�Jcst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		}
					 
		    node4 = node4 -> GetNext();
		   }
		    node3 = node3 -> GetNext();
	       }
           //�`�@M*N�����


//[15.5] ��V*u
/*	 
node2  = this-> GetNode_T();                    //��node����}����TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();            //j
	 
		node = this-> GetNode_B();              //��node����}����BS
			while(node != NULL){
			n = node->GetValue()->GetNo();      //i
				 
		node4 = this-> GetNode_T();             //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();    //j'
		    	
			if(n2!=n4)  // j!=j'
			{    	                
				for(int j=0 ; j<par ; j++)
				{
					Coes.add(1);
			        Vars.add(vars[ 2*M*N + 3*N + M ]);        //�s�Wu���ܼ�
			
			        Coes.add(-1);
			        Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);    //�s�W��V*u���ܼ�

			        Coes.add(Big);                                                                     //Big M
			        Vars.add(bilvars[ 2*N*M*N + 2*N*M*N*par + n2*M*N*par + n*N*par + n4*par + j ]);             //�s�W

			        sprintf(cst_name, "Mix_M_��V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);                      //��Mix_M_��V_%d_%d_%d_%d�g�Jcst_name
 			        IloConstraint cst;
			        cst = IloScalProd(Coes,Vars) <= Big;
			        cst.setName(cst_name);
			        model.add(cst); 
		   	        Coes.clear();
			        Vars.clear();	

		       }

		   }
					 
		     node4 = node4 -> GetNext();	
		    }
			 
		     node = node -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N�����

 
//[15.6] ��K*u
		 
node2  = this-> GetNode_T();                    //��node����}����TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();            //j
 		 
		node3 = this-> GetNode_B();             //��node����}����BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();    //i'
				 
		node4 = this-> GetNode_T();             //��node����}����TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();    //j'
		    	
		if(n2!=n4)    // j!=j' 
		{
			for(int j=0 ; j<par ; j++)
		    {
				Coes.add(1);
			    Vars.add(vars[ 2*M*N + 3*N + M ]);        //�s�Wu���ܼ�
			
			    Coes.add(-1);
			    Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);   //�s�W��K*u���ܼ�

			    Coes.add(Big);                                                                     //Big M
			    Vars.add(bilvars[ 2*N*M*N + 3*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);     //�s�W

			    sprintf(cst_name, "Mix_M_��K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                     //��Mix_M_��K_%d_%d_%d_%d�g�Jcst_name
 			    IloConstraint cst;
			    cst = IloScalProd(Coes,Vars) <= Big;
			    cst.setName(cst_name);
			    model.add(cst); 
		   	    Coes.clear();
			    Vars.clear();	

			}

		}		 
		     node4 = node4 -> GetNext();	
		    }
			 
		     node3 = node3 -> GetNext();		
		    }
			 
			 node2 = node2 -> GetNext();	
		    }
			//�`�@N*M*N�����
*/		   
 
//--------------------------------------------------------------------------------
 


IloCplex cplex(model);
//cplex.setParam(cplex.EpAGap,0.05);
cplex.solve();
//cplex.setParam(cplex.TiLim,10000); //����]���ɶ�
cplex.exportModel("BS10.lp");
cout<<"EE.lp"<<endl;

cout << "object finish." << endl; 
cout << "Solution status = " << cplex.getStatus() << endl;
cout << "Solution value  = " << cplex.getObjValue() << endl;  
	  
//---------------------------------------------------------------------------------------
	 	
		
		    char filename[250];
			FILE *fp;
			sprintf (filename,"BS10.txt");
			fp=fopen(filename,"w");
			double TT=0.0;
			double p_tol = 0.0;
 

//��ܭ�------------------------------------------------------------------------------------------------------
//Bin y	

		for(int i=0; i< M; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				  if(TT != 0)
					fprintf(fp,"\n[%s]:\t% .3f\n",binvars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//Bin x	
			
			for(int i=M; i< M + M*N ; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );
								 

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//Bin Zj
        	
			for(int i=M + M*N; i< M + M*N + N; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );
								 

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//Bin �fv
	    	
			for(int i=  M + M*N + N; i< M + M*N + N + M*N*par; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			
		fprintf(fp,"\n\n");

//Bin �fk      
		
			for(int i=  M + M*N + N + M*N*par; i< M + M*N + N + 2*M*N*par; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			
		fprintf(fp,"\n\n");
		fprintf(fp,"---------------------------------------------------------------\n");
		
//----------------------------------------------------------------------------------------------
//u
			
		for(int i= 2*M*N + 3*N + M; i< 1 + 2*M*N + 3*N + M; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];
				fprintf(fp,"\n[%s]:\t%.5f\n[1/%s]:\t%.5f\n",vars[i].getName(),TT,vars[i].getName() ,1/TT);
				p_tol = 1/TT;

				vals.clear();
				Vars.clear();

			}
			fprintf(fp,"\n\n");

//x*u	

		    fprintf(fp,"[x_u]\n");	
			for(int i= M*N + 2*N; i< M*N + 2*N + M*N; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .5f\n",vars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n"); 
		
//y*u

		    fprintf(fp,"[y_u]\n");	
			for(int i= 2*M*N + 2*N; i< 2*M*N + 2*N + M; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .5f\n",vars[i].getName(),TT  );


				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//�fv*u

		    fprintf(fp,"[�fv_u]\n");	
			for(int i= 1 + 2*M*N + 3*N + M; i< 1 + 2*M*N + 3*N + M + M*N*par ; i++){

				Vars.add( vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .5f\n",vars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//�fk*u

            fprintf(fp,"[�fk_u]\n");	
			for(int i= 1 + 2*M*N + 3*N + M + M*N*par; i< 1 + 2*M*N + 3*N + M + 2*M*N*par; i++){

				Vars.add( vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .5f\n",vars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//--------------------------------------------------------------------------	
//Power
			fprintf(fp,"[Power]\n");	
			for(int i=0; i< M*N; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
				//fprintf(fp,"[%s]:\t% .5f\n",vars[i].getName(),TT  );
				 fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",vars[i].getName() ,TT, p_tol*TT);
	


				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//----------------------------------------------------------------------------------------------
//SINR

		 	fprintf(fp,"[SINR]\n");	
			for(int i= M*N; i< M*N + N; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
				 fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",vars[i].getName() ,TT, p_tol*TT);

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

//----------------------------------------------------------------------------------------------
//H

			fprintf(fp,"[H = ln(1+SINR) ]\n");	
			for(int i= M*N +N; i < M*N + 2*N; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
				 fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",vars[i].getName() ,TT, p_tol*TT);

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");
 
//----------------------------------------------------------------------------------------------------------

 	fprintf(fp,"[ V ]\n");
		 	
	node2  = this-> GetNode_T();
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();
		
 		 
		node = this-> GetNode_B();
			while(node != NULL){
			n = node->GetValue()->GetNo();
		

		node4 = this-> GetNode_T();
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();
			 
					 
			if(n2!=n4){
			
				Vars.add(bilvars[ n2*M*N + n*N + n4 ]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
			 
                 fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",bilvars[ n2*M*N + n*N + n4 ].getName() ,TT, p_tol*TT);
	
				vals.clear();
				Vars.clear();
	
 	             }
			
			
			node4  = node4  ->GetNext();	
		    }
			 
		   node = node -> GetNext();		
		 }
		 
		node2 = node2 -> GetNext();	
	 }
			 
	 

fprintf(fp,"\n\n");

//-------------------------------------------------------------------------
fprintf(fp,"[ K ]\n");	
 
	node2  = this-> GetNode_T();
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();
		
 		 
		node3 = this-> GetNode_B();
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();
		

		node4 = this-> GetNode_T();
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();
			 
					 
			if(n2!=n4){
			
				Vars.add(bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
			   
				fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",bilvars[ N*M*N + n2*M*N + n3*N + n4 ].getName() ,TT, p_tol*TT);
				
				 vals.clear();
				 Vars.clear();
	
 	             }
					
			node4  = node4  ->GetNext();	
		    } 
			 
		   node3 = node3 -> GetNext();		
		 }
		 
		node2 = node2 -> GetNext();	
	 }
			 
	

	 fprintf(fp,"\n\n");

//-----------------------------------------------------------------------------------------------------------------------


fprintf(fp,"[��V]\n");	
	 		 
	node2  = this-> GetNode_T();
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();
	 		 
	node = this-> GetNode_B();
		while(node != NULL){
		n = node->GetValue()->GetNo();
		 

	node4 = this-> GetNode_T();
		while(node4 != NULL){
		n4 = node4->GetValue()->GetNo();
		 
		if(n2!=n4){ 
     		for(int j=0 ; j<par ; j++){
			

				Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
			    //fprintf(fp,"[%s]:\t %.5f\n",bilvars [2*M*N*M*N + n*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ].getName(),TT  );
					fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ].getName() ,TT, p_tol*TT);

				vals.clear();
				Vars.clear();
			 
			}}
			   node4  = node4  ->GetNext();	
		    }
			 
 
		     node = node -> GetNext();		
		    }
			 
 
			 node2 = node2 -> GetNext();	
		    }
			 	 
			 

	fprintf(fp,"\n\n");
//-----------------------------------------------------------------------------------------------------------
fprintf(fp,"[��K]\n");		 
			 
	node2  = this-> GetNode_T();
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();
	 		 
	node3 = this-> GetNode_B();
		while(node3 != NULL){
		n3 = node3->GetValue()->GetNo();
		 

	node4 = this-> GetNode_T();
		while(node4 != NULL){
		n4 = node4->GetValue()->GetNo();

		if(n2!=n4){ 
		for(int j=0 ; j<par ; j++){

				Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
			   // fprintf(fp,"[%s]:\t %.5f\n",bilvars [2*M*N*M*N + M*N*M*N*par + n*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ].getName(),TT  );
				fprintf(fp,"[%s]:\t%.5f\t\t%.5f\n",bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ].getName() ,TT, p_tol*TT);


				vals.clear();
				Vars.clear();
		}}
		  
			   node4  = node4  ->GetNext();	
		    }
			 
 
		     node3 = node3 -> GetNext();		
		    }
			 
 
			 node2 = node2 -> GetNext();	
		    }
			 


     	/*	 fprintf(fp,"[y]\n");	
			for(int i=0; i< M; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				 TT = vals[0];

				  if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

	//-----------------------------------------------------------------------------------------------------------		
			fprintf(fp,"[x]\n");	
			for(int i=M; i< M + M*N; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );
								 

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");


//--------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------------		
			fprintf(fp,"[z]\n");	
			for(int i=M + M*N; i< M +N+ M*N; i++){

				Vars.add(binvars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .3f\n",binvars[i].getName(),TT  );
								 

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n");

				fprintf(fp,"[P_D]\n");	
			for(int i=0; i< M*N; i++){

				Vars.add(vars[i]);
				cplex.getValues(vals, Vars);
				TT = vals[0];

				 if(TT != 0)
					fprintf(fp,"[%s]:\t% .7f\n",vars[i].getName(),TT  );

				vals.clear();
				Vars.clear();
			}
			fprintf(fp,"\n\n"); */





		 


		fprintf(fp,"\n\n");
		fprintf(fp,"---------------------------------------------------------------\n");

 
	fprintf(fp,"\n[ObjValue]:\t%.3f(%.3f)\n",cplex.getObjValue(),p_tol*cplex.getObjValue());
 
	
	
//----------------------------------------------------------------------------------------------
		

			fprintf(fp,"\n\nTime:%2f\n",cplex.getTime());
			fprintf(fp,"epagap:%2f\n",cplex.getParam(cplex.EpAGap));
			fprintf(fp,"epgap:%2f\n",cplex.getParam(cplex.EpGap));

			fclose(fp);	 

          
     // cplex.writeSolution("MILFP.sol");       

//---------------------------------------------
	
	   cplex.clear();
		vals.clear();
		Vars.clear();
		Coes.clear();
		vars.clear();
		model.end();


//--------------------------------------------------------------------------------
 
	 		
	 }
	catch (IloException& e){}//try 
  
	return;

 }//func

