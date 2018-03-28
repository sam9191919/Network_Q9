#pragma once 
#include "network.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <string>

#define mu	1.0                                                     //服務率
#define SF   1.0                                                    //Spreading factor，固定
#define Pn  (double) 4.11*pow(10.0,(-14))							//熱雜訊(mW)，注意單位，其中	Pn=kT k為Boltzmann常數	T為絕對溫度
#define Keq 2.75*pow(10.0,(15))										//Gain公式裡會用到
#define Alpha (float) 1.0											//正交因子，固定
#define W  10.0														//頻寬 , 如果這邊改了那注意熱雜訊Pn也要改
#define Ps 200.0													//基地台靜態功率W
#define bsmaxpower 40.0												//基地台所能提供功率總和上限W
#define k    1.0                                                    //使用者單位傳輸量所需支付的費用(dollars/bit)      //#define k  0.2
#define c    1.0                                                    //單位能量電費(dollars/joule)          //#define c	13.0/1000
#define pi   5.0                                                    //視利潤要為多少所設的變數
#define Big  (double) 20.0											//Big M
#define dmax 1500.0													//最大距離
#define beta  1.0													//目標式	

#define ga 6.0														//1為均勻切段，3~10 為非均勻切段(越大代表切越細)

#define par (int) 4													//bilinear term切段(切越多越複雜)##
#define nsum 4.0

//常用參數修改
#define PUB  600.0													//Power最大(mW)
#define SINR   100.0												//SINR最大
#define Rate 0.8													//datarate



using namespace std; 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
 

void NETWORK :: AddNode_T(NODE *node){                                           //加入新的TP
	if(Node_T==NULL){                                                            //如果原本的是空的，就是第一個
		Node_T = new LIST<NODE*>(node);                                          //new運算子會配置一個data type所需要的空間，並傳回該空間的位址；使用new運算子動態配置的空間，在整個程式結束前並不會自動歸還給記憶體，必須使用delete將這個空間還給記憶
		Node_T->GetValue()->SetNo(N++);                                          //填值後，序號++
		//cout<<"n="<< Node_T->GetValue()-> GetNo() <<endl;
		 
	}  else{                                                                     //如果不是空的
		Node_T->Tail()->SetNext(new LIST<NODE*>(Node_T->Tail(),node,NULL));      //從尾巴加入，依照(Node_T->Tail(),node,NULL)的順序
		Node_T->Tail()->GetValue()->SetNo(N++);                                  //填值後，序號++
	} 
}

void NETWORK :: AddNode_B(NODE *node){                                           //語法與TP相同 
	if(Node_B==NULL){
		Node_B = new LIST<NODE*>(node);
		Node_B->GetValue()->SetNo(M++);
		 
	} else{
		Node_B->Tail()->SetNext(new LIST<NODE*>(Node_B->Tail(),node,NULL));
		Node_B->Tail()->GetValue()->SetNo(M++);
	} 
} 




//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
void NETWORK :: AddUBLB(){                         //上、下限函式

int M,N,n,n2,n3,n4;                                //基地台數目M，使用者數目N

	LIST<NODE*> *node,*node2;                      //i,j
	LIST<NODE*> *node3,*node4;                     //i',j'
	

	M = this->GetNumberOfNode_B();                 //取得BS數量	 
	N = this->GetNumberOfNode_T();                 //取得TP數量

		//cout<<"UBLBHERE"<<endl;
		//system("pause");


//----------------------------------------------------------------------初始化-------------------------------------------------------------------------------------------------------//


	this->BinVar_LB = new int[ M + M*N + N + 2*M*N*par ];			 //Binary Variable，M為yi，M*N為xi,j，M*N*par為λv,λk i',j',np，N是zj
	this->BinVar_UB = new int[ M + M*N + N + 2*M*N*par ];
	for(int i=0;	  i< M + M*N + N + 2*M*N*par  ; i++){            //用迴圈初始化所有記憶體位址
		BinVar_LB[i] = 0;
		BinVar_UB[i] = 0;
	}

//------------------------------------------------------------------------------------------------

	this->Var_LB = new float[ 1 + 2*M*N + 3*N + M + 2*M*N*par ] ;         //Power,SINR,H=(1+SINR),只有Downlink，創建動態記憶體空間(有改)
	this->Var_UB = new float[ 1 + 2*M*N + 3*N + M + 2*M*N*par ] ;

	for(int i=0;   i< 1 + 2*M*N + 3*N + M + 2*M*N*par ; i++){             //用迴圈初始化所有記憶體位址
		Var_LB[i] = 0; 
		Var_UB[i] = 0;
	}

//-------------------------------------------------------------------------------------------------

	this->BilVar_LB = new float[  2*N*M*N + 2*N*M*N*par ]; //Biliear Tearm，2*N*M*N為K & V，；2*N*M*N*par為ΔK & ΔV，par指的是第n段
	this->BilVar_UB = new float[  2*N*M*N + 2*N*M*N*par ];

	for(int i=0;      i< 2*N*M*N + 2*N*M*N*par ; i++){     //用迴圈初始化所有記憶體位址
		 BilVar_LB[i] = 0; 
		 BilVar_UB[i] = 0; 
	}
 
//------------------------------------------------------------變數上下限------------------------------------------------//
 //u
 
     cout<<"Connection = 0.8"<<endl;
     cout<<"Rate = 0.35M"<<endl;
     cout<<"Ps = 200W "<<endl;
    

	//Var_UB[ 1+ M + 4*M*N ] =1/(20*M*N + (Ps*M));				//u 分數計算的部分

	 Var_UB[ 2*M*N + 3*N + M  ] = 2.0;							//u的Upper Bound
    //cout<<"u="<<Var_UB[ 2*M*N + 3*N + M  ]<<endl;

     Var_LB[ 2*M*N + 3*N + M  ] = 0;							//u的Lower Bound


//------------------------------------------------------------------------Binary Term && Variable Term---------------------------------------------//

node = this-> GetNode_B();                                                           //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                                   //i

	node2  = this-> GetNode_T();                                                     //用node的位址找到該TP
	while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                                             //j
				
//Varible 宣告
						
		
		        Var_UB[ 	    n*N + n2  ] = 250 * Var_UB[ 2*M*N + 3*N + M  ];                           //Power
				Var_UB[               M*N + n2  ] = SINR  * Var_UB[ 2*M*N + 3*N + M  ];	                  //SINR*u 
				Var_UB[           M*N + N + n2  ] = log(1+Var_UB[ M*N + n2 ]);                            //H = (1+SINR)		
                Var_UB[   M*N + 2*N + n*N + n2  ] = Var_UB[ 2*M*N + 3*N + M ];                            //x*u = u
                Var_UB[       2*M*N + 2*N + n   ] = Var_UB[ 2*M*N + 3*N + M ];                            //y*u = u
				Var_UB[   2*M*N + 2*N + M + n2  ] = Var_UB[ 2*M*N + 3*N + M ];                            //Zj*u = u
			//	Var_UB[       2*M*N + 3*N + M   ] = Var_UB[ 2*M*N + 3*N + M ];                            //u





//Binary 宣告		
		        BinVar_UB [           n   ] = 1.0;                                            //yi
				BinVar_UB [ M + n*N + n2  ] = 1.0;                                            //xij
				BinVar_UB [ M + M*N + n2  ] = 1.0;                                            //Zj 			  



			node2  = node2  ->GetNext();	
		} 
		 node  = node  ->GetNext();	
	}



//------------------------------------------------------------------------Bilinear Term---------------------------------------------//

//( V & K )*U

	 		 
		node = this -> GetNode_B();																				//用node的位址找到該BS
	while ( node != NULL )
	{
		n = node -> GetValue() -> GetNo();																		//i
		node2 = this -> GetNode_T();																			//用node的位址找到該TP
		
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


			                                   BilVar_UB[     n2*M*N + n*N + n4 ] = 100 * SINR * Var_UB[ 2*M*N + 3*N + M ] ;                           //V對應到的記憶體位置 => Powerij'*SINRj*U
					                       
											   BilVar_UB[     N*M*N + n2*M*N + n3*N + n4 ] = 100 * SINR * Var_UB[ 2*M*N + 3*N + M ] ;                  //K對應到的記憶體位置 => Poweri'j'*SINRj*U	                               
				 
				 node4 = node4 -> GetNext();	
		      } 
			  node3 = node3 -> GetNext();		
		   } 
			node2 = node2 -> GetNext();	
		} 
		node  = node -> GetNext();	
	} 
//△V & △K


	node = this -> GetNode_B();																					//同上
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
					                           	 BilVar_UB[                 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ] = PUB * SINR * Var_UB[ 2*M*N + 3*N + M ] ;           //△V = Powerij'*SINRj*U的upper bound

												 BilVar_UB[     2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ] = PUB * SINR * Var_UB[ 2*M*N + 3*N + M ] ;         //△K = Poweri'j'*SINRj*U的upper bound	 
  
				                         
				                              
			    
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

//λv && λv*u
node = this-> GetNode_B();                                                            //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                                    //i

	node4  = this-> GetNode_T();                                                      //用node的位址找到該TP
	    while(node4 != NULL){
		n4 = node4->GetValue()->GetNo(); 	                                          //j'	                                          

             for(int j=0; j<par; j++)
			 {
                                 
				   BinVar_UB[       M + M*N + N + n*N*par + n4*par + j  ] = 1 ;                                         //將在相同的n及n4下切段的第j段設為λv
				   Var_UB[  1 + 2*M*N + 3*N + M + n*N*par + n4*par + j  ] = Var_UB[ 2*M*N + 3*N + M ];                  //λv*u = u			   
			 } 
		  
		node4 = node4 -> GetNext();	                                                  //指向下一個TP
	   }  
	node = node -> GetNext();	                                                      //指向下一個BS
   } 

//λk  && λk*u
 node3 = this-> GetNode_B();                                                           //用node的位址找到該BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();                                                   //i'

	node4 = this-> GetNode_T();                                                        //用node的位址找到該TP
	    while(node4 != NULL){
		n4 = node4->GetValue()->GetNo(); 	                                           //j'	                                           

             for(int j=0; j<par; j++)
			 {
                                 		   
				   BinVar_UB[       M + M*N + N + M*N*par + n3*N*par + n4*par + j  ] = 1 ;                               //將在相同的n3及n4下切段的第j段設為λk
				   Var_UB[  1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j  ] = Var_UB[ 2*M*N + 3*N + M ];        //λk*u = u
			 } 
		        
     	node4 = node4 -> GetNext();	                                                   //指向下一個TP
	   }  
	node3 = node3 -> GetNext();	                                                       //指向下一個BS
   }

    return;
	system("pause");

	}

//----------------------------Add Variables 位置---------------------------------------------------------//

void NETWORK ::func(){ //NETWORK *network

	//cout<<"Vars Location"<<endl;


    int M,N;

	M = this->GetNumberOfNode_B();	 
	N = this->GetNumberOfNode_T();
	
	ILOSTLBEGIN                                                          //cplex的標準函式庫
 
	
	IloEnv env;                                                          //建置env環境
	IloModel model(env);                                                 //env環境裡建立model模型
	IloNumVarArray vars(env), Vars(env), binvars(env), bilvars(env);     //宣告陣列
	IloNumArray Coes(env), vals(env);                                    //宣告數值
	
	

	LIST<NODE*> *node, *node2;                                          //i,j
	LIST<NODE*> *node3,*node4;                                          //i',j'
	LIST<NETWORK*> *network;											//這裡重複宣告是為了將位址導入cplex

	int n,n2,n3,n4,i,j;	

 
  try{                                                                   //錯誤偵測語法

		//if(BNB_DDisplay){
		//	cout << "Adding vars."<<endl;

	
 char var_name[60];                                                      //設一個大小為60的陣列,存變數名字
// char cst_name[60];                                                    //設一個大小為60的陣列,存常數名字
		 

// Add Power_D
node  = this-> GetNode_B();                                              //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                       //i

    node2  = this-> GetNode_T();                                         //用node2的位址找到該TP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo();                             //j
				
				char var_name[60];                                                                             //設一個大小為60的陣列 ##
				sprintf(var_name,"P_%d_%d",n+1,n2+1);                                                          //把P_D_%d_%d寫入var_name，並印出該參數i,j
			 	vars.add(IloNumVar(env, Var_LB[ n*N + n2 ], Var_UB[ n*N + n2 ], var_name));                    //把c加到cplex，power_D要在上下限之間

					  
			 node2 = node2 -> GetNext();	
		    } 
			 
 
		node  = node  -> GetNext();	
		}  


//Add SINR_D                                           
	
    node2  = this-> GetNode_T();                                           //與Power_D語法相同
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo(); 
			
				char var_name[60];                                                                              
				sprintf(var_name,"S_%d" ,n2+1);                                                                
			 	 vars.add(IloNumVar(env, Var_LB[ M*N + n2 ], Var_UB[ M*N + n2 ], var_name));
				
				 node2 = node2 -> GetNext();	
		        }
	 
			  

//Add H_D = 1+SINR

    node2  = this-> GetNode_T();                                           //與Power_D語法相同
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo(); 

			char var_name[60];
			sprintf(var_name,"H_%d",n2+1);
		 	vars.add(IloNumVar(env, Var_LB[ M*N + N + n2 ], Var_UB[ M*N + N + n2 ], var_name));
 
			     node2 = node2 -> GetNext();	
		        }
			 

//Add x*u

node  = this-> GetNode_B();																			//用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();																	//i

    node2  = this-> GetNode_T();																	//用node的位址找到該TP
		while (node2 != NULL){  
	    n2 = node2->GetValue()->GetNo();															//j

			char var_name[60];
			sprintf(var_name,"x_%d_%d",n+1,n2+1);
		 	vars.add(IloNumVar(env, Var_LB[ M*N + 2*N + n*N + n2 ], Var_UB[ M*N + 2*N + n*N + n2 ], var_name));   //x*u的記憶體位置
			 
	
			     node2 = node2 -> GetNext();	
		        } 
			 
	 
			 node  = node  -> GetNext();	
		    }  

//Add y*u

node  = this-> GetNode_B();																			//用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();																	//i

			char var_name[60];
			sprintf(var_name,"y_%d",n+1 );
		 	vars.add(IloNumVar(env, Var_LB[ 2*M*N + 2*N + n ], Var_UB[ 2*M*N + 2*N + n ], var_name));                  //y*u的記憶體位置
 
			 node  = node  -> GetNext();	
		    }  
 
//Add Zj*u

node2  = this-> GetNode_T();																		//用node的位址找到該BS
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();																//j

			char var_name[60];
			sprintf(var_name,"Zj_%d",n2+1 );
		 	vars.add(IloNumVar(env, Var_LB[ 2*M*N + 2*N + M + n2 ], Var_UB[ 2*M*N + 2*N + M + n2 ], var_name));        //Zj*u的記憶體位置
 
			 node2  = node2  -> GetNext();	
		    } 

	//Add u
	
	//char var_name[60];
	sprintf(var_name,"u");
	vars.add(IloNumVar(env, Var_LB[ 2*M*N + 3*N + M ], Var_UB[ 2*M*N + 3*N + M ], var_name));                          //u的記憶體位置


//Add λv*u

node  = this-> GetNode_B();                                                                         //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();																	//i

    node4  = this-> GetNode_T();																	//用node的位址找到該TP
		while (node4 != NULL){  
			n4 = node4->GetValue()->GetNo();														//j'

				for(int j=0 ; j<par ; j++)
				{
					char var_name[60];
					sprintf(var_name,"λv_%d_%d_%d",n+1,n4+1,j+1);
		 			vars.add(IloNumVar(env, Var_LB[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ], Var_UB[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ], var_name));                       //λv*u的記憶體位置
			    }
	
			  node4 = node4 -> GetNext();	
		     }  
			 
 
		node  = node  -> GetNext();	
	   }  

//Add λk*u

node3  = this-> GetNode_B();                                                                          //用node的位址找到該BS
	while(node3 != NULL){                                                                             
	n3 = node3->GetValue()->GetNo();                                                                  //i'

    node4  = this-> GetNode_T();                                                                      //用node的位址找到該TP
		while (node4 != NULL){  
			n4 = node4->GetValue()->GetNo();                                                          //j'

				for(int j=0 ; j<par ; j++)
				{
					char var_name[60];
					sprintf(var_name,"λk_%d_%d_%d",n3+1,n4+1,j+1);
		 			vars.add(IloNumVar(env, Var_LB[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ], Var_UB[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ], var_name));    //λk*u的記憶體位置
			    }
	
			  node4 = node4 -> GetNext();	
		     }  
			 
 
	      node3  = node3  -> GetNext();	
		 }




//------------------------------------------------------Add Binary 位置----------------------------------------------------------------------//

//Add y
node  = this-> GetNode_B();                                                                           //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                                                    //i

	      		char var_name[60];                                                                    //宣告變數名稱的記憶體 ##
				sprintf(var_name,"Bin_y_%d",n+1);                                                     //把Bin_y_%d寫入var_name，並印出i,j
			 	binvars.add(IloNumVar(env, BinVar_LB[ n ], BinVar_UB[ n ],ILOINT , var_name));        //將y加入cplex環境中
				//cout <<"y:BinVar_LB = "<<i<< endl;
		   	  node  = node  -> GetNext();	
		    } 
 
//Add x
node  = this-> GetNode_B();                                                                                       //用node的位址找到該BS
	while(node != NULL){                                                                         
	n = node->GetValue()->GetNo();                                                                                //i
 
    node2  = this-> GetNode_T();                                                                                  //用node的位址找到該TP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo();                                                                      //j

			char var_name[60];                                                                                    //宣告變數名稱的記憶體 ##
			sprintf(var_name,"Bin_x_%d_%d",n+1,n2+1);                                                             //把Bin_x_%d_%d寫入var_name，並印出i,j
		 	binvars.add(IloNumVar(env, BinVar_LB[ M + n*N + n2 ], BinVar_UB[ M + n*N + n2 ],ILOINT , var_name));  //將x加入cplex環境中
 
			 
			   node2 = node2 -> GetNext();	
		      }
			 
		   node  = node  -> GetNext();
	      }

	
//Add Zj
node2  = this-> GetNode_T();                                                                                                //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                                                                        //j

	      		char var_name[60];                                                                                          //宣告變數名稱的記憶體 ##
				sprintf(var_name,"Bin_Zj_%d",n2+1);                                                                         //把Bin_Zj_%d寫入var_name，並印出j
			 	binvars.add(IloNumVar(env, BinVar_LB[ M + M*N + n2 ], BinVar_UB[ M + M*N + n2 ],ILOINT , var_name));        //將Zj加入cplex環境中
	
		   	  node2  = node2  -> GetNext();	
		    } 
				
//Addλv  	  
node  = this-> GetNode_B();                                                                                                                                        
	while(node != NULL){ 
	n = node->GetValue()->GetNo();                                                                                
 
    node4  = this-> GetNode_T();                                                                                   
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                                                       

						for(int j=0 ; j<par ; j++)
						{

							char var_name[60];
							sprintf(var_name,"Bin_λv_%d_%d_%d ",n+1,n4+1,j+1);                                                                                                        //把Bin_λv_%d_%d_%d寫入var_name，並印出i,j,par
		 					binvars.add(IloNumVar(env, BinVar_LB[ M + M*N + N + n*N*par + n4*par + j ], BinVar_UB[ M + M*N + N + n*N*par + n4*par + j ],ILOINT , var_name));           //將λv加入cplex環境中
						   
						}
 				  node4 = node4-> GetNext();
				 }   
	        
			node = node -> GetNext();
	       }

//Addλk  	  
node3  = this-> GetNode_B();                                                                                                                                        
	while(node3 != NULL){ 
	n3 = node3->GetValue()->GetNo();                                                                                
 
    node4  = this-> GetNode_T();                                                                                   
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                                                       

						for(int j=0 ; j<par ; j++)
						{

							char var_name[60];
							sprintf(var_name,"Bin_λk_%d_%d_%d ",n3+1,n4+1,j+1);                                                                                                                        //把Bin_λk_%d_%d_%d寫入var_name，並印出i,j,par
		 					binvars.add(IloNumVar(env, BinVar_LB[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ], BinVar_UB[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ],ILOINT , var_name));      //將λk加入cplex環境中
						   
						}
 				  node4 = node4-> GetNext();
				 }   
	        
			node3 = node3 -> GetNext();
	       }

 
//-----------------------------------------------------Add Bilinear Term 位置------------------------------------------------------//
 
//Add V=Powerij'*SINRj*U		

	node2  = this-> GetNode_T();                                                                                                                                //用node的位址找到該TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                                                                                                                        //j
			//cout<< "n2=" << node2->GetValue()->GetNo() <<endl;
 		   node = this-> GetNode_B();                                                                                                                           //用node的位址找到該BS
			  while(node != NULL){
			  n = node->GetValue()->GetNo();                                                                                                                    //i
				//cout<< "n3=" << node3->GetValue()->GetNo() <<endl;
                 node4 = this-> GetNode_T();                                                                                                                    //用node的位址找到該TP
		              while(node4 != NULL){
			          n4 = node4->GetValue()->GetNo();                                                                                                          //j'
			  
			                      char var_name[60];                                                                                                            //宣告變數名稱的記憶體 ##
			                      sprintf(var_name,"V_%d_%d_%d",n2+1,n+1,n4+1);                                                                                 //把V_%d_%d_%d寫入var_name，並印出j,i,j'
		 	                      bilvars.add(IloNumVar(env, BilVar_LB[ n2*M*N + n*N + n4 ], BilVar_UB[ n2*M*N + n*N + n4 ], var_name));                        //將V加入cplex環境中
			
			   
			     node4  = node4  ->GetNext();	
		        }
			 
		    node = node -> GetNext();		
		   }
		
		node2 = node2 -> GetNext();	
	   }
			 


//Add K=Poweri'j'*SINRj*U

	node2  = this-> GetNode_T();                                                                                                                                        //用node的位址找到該TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                                                                                                                                //j
			//cout<< "n2=" << node2->GetValue()->GetNo() <<endl;
 		   node3 = this-> GetNode_B();                                                                                                                                  //用node的位址找到該BS
			   while(node3 != NULL){
			   n3 = node3->GetValue()->GetNo();                                                                                                                         //i'
				//cout<< "n3=" << node3->GetValue()->GetNo() <<endl;
                  node4 = this-> GetNode_T();                                                                                                                           //用node的位址找到該TP
		              while(node4 != NULL){
			          n4 = node4->GetValue()->GetNo();                                                                                                                  //j'

                                char var_name[60];                                                                                                                      //宣告變數名稱的記憶體 ##
 			                    sprintf(var_name,"K_%d_%d_%d",n2+1,n3+1,n4+1 );                                                                                         //把K_%d_%d_%d寫入var_name，並印出j,i',j'
	 		                    bilvars.add(IloNumVar(env,BilVar_LB[ N*M*N + n2*M*N + n3*N + n4 ], BilVar_UB[ N*M*N + n2*M*N + n3*N + n4 ], var_name));                 //將K加入cplex環境中
			                 
								

			  node4  = node4  ->GetNext();	
		     }     
			 
		   node3 = node3 -> GetNext();		
		  }
		
		node2 = node2 -> GetNext();	
	   } 
			 


 //Add △V
			 
	node2  = this-> GetNode_T();                                 //用node的位址找到該TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                         //j
	 		 
	node = this-> GetNode_B();                                   //用node的位址找到該BS
		while(node != NULL){
		n = node->GetValue()->GetNo();                           //i
		 
	node4 = this-> GetNode_T();                                  //用node的位址找到該TP
		while(node4 != NULL){
		n4 = node4->GetValue()->GetNo();                         //j'
		 
		for(int j=0 ; j<par ; j++)
		{

 				char var_name[60];
				sprintf(var_name,"△V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);                                                                                                         //△V名稱一樣注意 ##
	 			bilvars.add(IloNumVar(env, BilVar_LB[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ], BilVar_UB[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ], var_name));     //將△V加入cplex環境中
			    
				//cout << "Sijn: BilVar_LB["<<i<<"*"<<M<<"*"<<N<<"*"<<par<<"+"<< i<<"*"<<N<<"*"<<par<<"+"<<n2<<"*"<<par<<"+"<<j<<"]="<< 4*M*N*M*N +n*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j<<endl;
	    }
	
			 node4 = node4 -> GetNext();	
		    }
			 
 
		     node = node -> GetNext();		
		    }
			 
 
			 node2 = node2 -> GetNext();	
		    }
			 
		 
			
//Add △K
			 
	node2  = this-> GetNode_T();                                  //用node的位址找到該TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();                          //j
	 		 
	node3 = this-> GetNode_B();                                   //用node的位址找到該BS
		while(node3 != NULL){
		n3 = node3->GetValue()->GetNo();                          //i'
		 

	node4 = this-> GetNode_T();                                   //用node的位址找到該TP
		while(node4 != NULL){
		n4 = node4->GetValue()->GetNo();                          //j'
		 
		for(int j=0 ; j<par ; j++)
		{

 				char var_name[60];
				sprintf(var_name,"△K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                                                                                                                                 //△K名稱一樣注意 ##
	 			bilvars.add(IloNumVar(env, BilVar_LB[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ], BilVar_UB[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ], var_name));    //將△K加入cplex環境中
			    
				//cout << "Sijn: BilVar_LB ="<<4*M*N*M*N + N*M*M*N*par + n*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j  <<endl;

		}

			 node4 = node4 -> GetNext();	
		    }
			 
		   
		     node3 = node3 -> GetNext();		
		    }
			 
 
			 node2 = node2 -> GetNext();	
		    }	 
			


			 
 
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^---變數宣告結束--------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 	//Object function ##############

 
	for(int i= 2*M*N + 2*N + M; i< 2*M*N + 3*N + M ; i++){
			Vars.add(vars[i]);             //Zj*u的記憶體位置
			Coes.add(Rate*k );          //R*b
		 	} 

 
		model.add(IloMaximize(env, IloScalProd(Coes,Vars)));

		Coes.clear();
		Vars.clear(); 

	  
//----------------------------------目標函式宣告結束-------------------------------	

//----------------------------------------[Ij][Ji]矩陣宣告------------------------------------------

	double d;
	double Ij[200][200];																//第j個TP可以連的所有BS
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

				Ij[n2][n] = n + 1;														//若距離在dmax內 則填入BS編號

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
	double Jicounter[200];																			//第i個BS可以連的所有TP

	node = this->GetNode_B();
	while (node != NULL) {
		n = node->GetValue()->GetNo();
		//cout<< "n="<<node->GetValue()->GetNo()<<endl;

		node2 = this->GetNode_T();
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();
			//cout<< "n2="<<node2->GetValue()->GetNo()<<endl;

			d = sqrt(pow(node2->GetValue()->GetX1() - node->GetValue()->GetX2(), 2) + pow(node2->GetValue()->GetY1() - node->GetValue()->GetY2(), 2) + pow((node2->GetValue()->GetZ1()) - node->GetValue()->GetZ2(), 2));
			//TP與BS間距離

			if (d <= dmax)
			{

				Ji[n][n2] = n2 + 1;																	//若距離在dmax內 則填入TP編號

				Jicounter[n] = Jicounter[n] + 1;

				cout << "Ji" << n << n2 << "=" << Ji[n][n2] << endl;

				cout << "Jicounter" << n << "=" << Jicounter[n] << endl;


			}


			node2 = node2->GetNext();

		}
		node = node->GetNext();
	}

	//----------------------------------------------------------------------------------
	//儲存gain矩陣

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




			Hnorm = (137.4 + 35.2*(log10((d / 1000.0)))) / 10.0;//Hnorm=Lij/10   ;路徑損耗Lij   ;d單位為公里
			gain[n][n2] = 1 / (pow(10.0, Hnorm));//通道增益

												 /*  Hnorm = (40*(log10(d/1000.0))+30*(log10(fa))+49)/10.0;




												 gain[n][n2] = 1/ (pow(10,Hnorm)) ; */

			cout << "gain" << n << n2 << "=" << gain[n][n2] << endl;
			cout << "d" << n << n2 << "=" << d << endl;

			node2 = node2->GetNext();

		}
		node = node->GetNext();
	}

	//---------------------------------------限制式-------------------------------------------

	//[A]分數規劃令新變數u後多出的限制式
	node = this->GetNode_B();                 //用node的位址找到該BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();          //i

		node2 = this->GetNode_T();                //用node的位址找到該TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();        //j

			if (n == ((Ij[n2][n]) - 1))
			{
				if (n2 == ((Ji[n][n2]) - 1))
				{
					Coes.add(beta);
					Vars.add(vars[n*N + n2]);   //新增Power*u的變數  
				}                                   //Coes與Vars陣列為暫存係數與變數 使用完清空 
			}
			node2 = node2->GetNext();
		}                                           //beta*ΣΣPij*u

		if (n == ((Ij[n2][n]) - 1))
		{
			Coes.add(Ps);                                        //Ps
			Vars.add(vars[2 * M*N + 2 * N + n]);                  //y*u
		}


		node = node->GetNext();
	}                                           //ΣPs*yi*u

	char cst_name[60];
	sprintf(cst_name, "u");//把u寫入cst_name
	IloConstraint cst;//創建限制式cst
	cst = IloScalProd(Coes, Vars) == 1;//限制式A
	cst.setName(cst_name);//設立限制式名稱 u
	model.add(cst); //將限制式加入模型中
	Coes.clear();//清空Coes陣列
	Vars.clear();//清空Vars陣列


				 //[1]BS沒有開啟，則TP無法與該BS進行連結
	node = this->GetNode_B();                          //用node的位址找到該BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();                    //i

		node2 = this->GetNode_T();                      //用node的位址找到該TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();              //j


			Coes.add(1);
			Vars.add(vars[M*N + 2 * N + n * N + n2]);   //新增x*u的變數
			Coes.add(-1);
			Vars.add(vars[2 * M*N + 2 * N + n]);        //新增y*u的變數



			sprintf(cst_name, "BS_%d_%d", n + 1, n2 + 1); //把BS_%d_%d寫入cst_name
			IloConstraint cst;                    //創建限制式cst
			cst = IloScalProd(Coes, Vars) <= 0;     //限制式1 xij*u - yi*u <= 0
			cst.setName(cst_name);                //設立限制式名稱BS_%d_&d
			model.add(cst);                       //將限制式加入模型中
			Coes.clear();                         //清空Coes陣列
			Vars.clear();                         //清空Vars陣列

			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}
	//共M*N條限制式


	//[2]我們假設基地台的最大功率上限為Pmax同時，當TP沒有連結到某個BS時，該BS也不會分配下行傳輸功率給該TP。

	node = this->GetNode_B();                           //用node的位址找到該BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();                    //i

		node2 = this->GetNode_T();                      //用node的位址找到該TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();              //j

			Coes.add(-1.0);
			Vars.add(vars[n*N + n2]);           //新增Power*u的變數

			Coes.add(Var_UB[n*N + n2]);
			Vars.add(vars[M*N + 2 * N + n * N + n2]);  //新增x*u的變數


			char cst_name[60];
			sprintf(cst_name, "Power_%d.%d", n + 1, n2 + 1); //把Power_%d.%d寫入cst_name
			IloConstraint cst;                    //創建限制式cst
			cst = IloScalProd(Coes, Vars) >= 0;    //生成一條限制式
			cst.setName(cst_name);                //設立限制式名稱Power_%d
			model.add(cst);                       //將限制式加入模型中
			Coes.clear();                         //清空Coes陣列
			Vars.clear();                         //清空Vars陣列

			node2 = node2->GetNext();
		}
		node = node->GetNext();
	}
	//共M*N條限制式


	//[3]BS也必須是在開啟的情況下，才可以服務所連結的TP。且BS分配給所有TP的功率總合也不能超過Pmax

	node = this->GetNode_B();                           //用node的位址找到該BS
	while (node != NULL) {
		n = node->GetValue()->GetNo();                     //i

		node2 = this->GetNode_T();                       //用node的位址找到該TP
		while (node2 != NULL) {
			n2 = node2->GetValue()->GetNo();               //j

			Coes.add(1);
			Vars.add(vars[M*N + 2 * N + n * N + n2]);      //新增Power*u的變數

			node2 = node2->GetNext();
		}
		Coes.add(-bsmaxpower);                           //BS的Max Power
		Vars.add(vars[2 * M*N + 2 * N + n]);           //新增y*u的變數  


		char cst_name[60];
		sprintf(cst_name, "SumPower_%d", n + 1); //把SumPower_%d寫入cst_name
		IloConstraint cst;                     //創建限制式cst
		cst = IloScalProd(Coes, Vars) <= 0;      //生成一條限制式
		cst.setName(cst_name);                 //設立限制式名稱SumPower_%d
		model.add(cst);                        //將限制式加入模型中
		Coes.clear();                          //清空Coes陣列
		Vars.clear();                          //清空Vars陣列

		node = node->GetNext();
	}
	//共M條限制式


	//[4.1]前
	//若Ej=0(未進入服務) 服務其基地台數目為0
	//若Ej=1則最少可連1個基地台 最多可以連Ij台基地台

	node2 = this->GetNode_T();                          //用node的位址找到該TP
	while (node2 != NULL) {
		n2 = node2->GetValue()->GetNo();              //j

		node = this->GetNode_B();                       //用node的位址找到該BS
		while (node != NULL) {
			n = node->GetValue()->GetNo();                //i
			if ((n == ((Ij[n2][n]) - 1)))
			{
				Coes.add(1.0);
				Vars.add(vars[M*N + 2 * N + n * N + n2]);   //新增x*u的變數
			}

			node = node->GetNext(); //Σxij*u
		}

		Coes.add(-1.0);
		Vars.add(vars[2 * M*N + 2 * N + M + n2]);     //新增Zj*u的變數


		char cst_name[60];
		sprintf(cst_name, "TP-1-0_%d", n2 + 1);   //把TP_%d寫入cst_name
		IloConstraint cst;                    //創建限制式cst
		cst = IloScalProd(Coes, Vars) >= 0;     //生成一條限制式 Σxij*u - zj*u >=0
		cst.setName(cst_name);                //設立限制式名稱TP_%d
		model.add(cst);                       //將限制式加入模型中
		Coes.clear();                         //清空Coes陣列
		Vars.clear();	                      //清空Vars陣列

		node2 = node2->GetNext();
	}


	//[4.2]後

	node2 = this->GetNode_T();                          //用node的位址找到該TP
	while (node2 != NULL) {
		n2 = node2->GetValue()->GetNo();               //j

		node = this->GetNode_B();                        //用node的位址找到該BS
		while (node != NULL) {
			n = node->GetValue()->GetNo();                 //i
			if ((n == ((Ij[n2][n]) - 1)))
			{
				Coes.add(1.0);
				Vars.add(vars[M*N + 2 * N + n * N + n2]);    //新增x*u的變數
			}

			node = node->GetNext(); //Σxij*u
		}

		Coes.add(-(Ijcounter[n2]));
		Vars.add(vars[2 * M*N + 2 * N + M + n2]);      //新增Zj*u的變數


		char cst_name[60];
		sprintf(cst_name, "TP-1-1_%d", n2 + 1);    //把TP_%d寫入cst_name
		IloConstraint cst;                     //創建限制式cst
		cst = IloScalProd(Coes, Vars) <= 0;      //生成一條限制式 Σxij*u - Ij*Zj*u <= 0
		cst.setName(cst_name);                 //設立限制式名稱TP_%d
		model.add(cst);                        //將限制式加入模型中
		Coes.clear();                          //清空Coes陣列
		Vars.clear();	                       //清空Vars陣列

		node2 = node2->GetNext();
	}											//共N條限制式
  


//[5]為了滿足用戶對系統的使用滿意度，在整個行動網路覆蓋範圍內，至少要服務 mu%以上的TP服務數量，其中，mu為已知常數，代表最低Service Rate 之要求。 

node2  = this-> GetNode_T();                        //用node的位址找到該TP
	  while (node2 != NULL){  
	  n2 = node2->GetValue()->GetNo();              //j
	       
				Coes.add(1);
				Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);    //新增Zj*u的變數					

			node2  = node2 ->GetNext();	
		   }	
		
			Coes.add(-(N*mu));
			Vars.add(vars[ 2*M*N + 3*N + M ]);      //新增u的變數

 
			char sct_name[60];
			sprintf(sct_name, "outtage ");          //把outtage 寫入sct_name
			IloConstraint sct;
			sct = IloScalProd(Coes,Vars) >= 0 ;     //生成一條限制式
			sct.setName(sct_name);                  //設立限制式名稱outtage 
		    model.add(sct);                         //將限制式加入模型中
			Coes.clear();                           //清空Coes陣列
			Vars.clear();                           //清空Vars陣列


//[6]Distance 其中dmax為BS(i)和TP(j)之間最大可連線距離，是個已知常數。 

			double Gain;

	node  = this-> GetNode_B();                            //用node的位址找到該BS
		while(node != NULL){
		n = node->GetValue()->GetNo();                      //i

		node2  = this-> GetNode_T();                        //用node的位址找到該TP
			while (node2 != NULL){  
				n2 = node2->GetValue()->GetNo();            //j
				
				d = sqrt(pow(node2->GetValue()->GetX1()-node ->GetValue()->GetX2(),2)+pow(node2->GetValue()->GetY1()-node ->GetValue()->GetY2(),2)+pow((node2->GetValue()->GetZ1())-node ->GetValue()->GetZ2(),2));				//距離為平方開根號	
							
				Coes.add(d);
				Vars.add(vars[ M*N + 2*N + n*N + n2 ]);      //x*u

				Coes.add(-dmax);
				Vars.add(vars[ 2*M*N + 3*N + M ]);           //新增u的變數

					char var_name[60];
					sprintf(cst_name, "Distance_%d_%d", n+1, n2+1);     //把Distance_%d.%d寫入cst_name
					IloConstraint cst;                                  //創建限制式cst
					cst = IloScalProd(Coes,Vars) <=0;                   //生成一條限制式
					cst.setName(cst_name);                              //設立限制式名稱Distance_%d.%d
					model.add(cst);                                     //將限制式加入模型中
					Coes.clear();                                       //清空Coes陣列
					Vars.clear();                                       //清空Vars陣列
			
				node2  = node2  ->GetNext();
                 	
				}
			node  = node  ->GetNext();
			} 
		 
 			//總共M*N條限制式

 //[7]   //d2 會不會也要綁  //改到這裡
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

								
					d = sqrt(pow(node2->GetValue()->GetX1()-node ->GetValue()->GetX2(),2)+pow(node2->GetValue()->GetY1()-node ->GetValue()->GetY2(),2)+pow((node2->GetValue()->GetZ1())-node ->GetValue()->GetZ2(),2));							//距離為平方開根號



						/*	Hnorm= pow (10.0, 3.25*(log10(d/1000.0)) );
							Gain = 1.0/ (Keq*Hnorm) ;*/

	
					/*    Hnorm = (40*(log10(d/1000.0))+30*(log10(2300.0))+49)/10.0;
							  Gain = 1/ (pow(10,Hnorm)) ; */


							Hnorm = (137.4 + 35.2*(log10((d/1000.0))))/10.0;				//Lij(路徑損耗)/10	
							Gain = 1.0/ (pow(10.0,Hnorm)) ;									//clannel gain

						/*	Hnorm = (143.503 + 38.35*(log10(d/1000.0)))/10.0;
							Gain = 1/ (pow(10,Hnorm)) ; */

		



						   d2 = sqrt(pow(node2->GetValue()->GetX1()-node3 ->GetValue()->GetX2(),2)+pow(node2->GetValue()->GetY1()-node3 ->GetValue()->GetY2(),2)+pow((node2->GetValue()->GetZ1())-node3 ->GetValue()->GetZ2(),2));				//距離為平方開根號



						/*		Hnorm2= pow (10.0, 3.25*(log10(d2/1000.0)) );
								Gain2 = 1.0/(Keq*Hnorm2) ;*/

			
						/*		Hnorm2 = (40*(log10(d2/1000.0))+30*(log10(2300.0))+49.0)/10.0;
			

			

								Gain2 = 1/ (pow(10,Hnorm2)) ;*/
								Hnorm2 = (137.4 + 35.2*(log10((d2/1000.0))))/10.0;			//同上
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
											 if ((n!=((Ij[n2][n])-1))&&(n4==(((Ji[n3][n4])-1)))){				//Ij是第j個TP可以連上的所有BS , Ji是第i個BS可以連上的所有TP	
																												//J不在I的覆蓋範圍內，J'在I的覆蓋範圍內
		        
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
										 if((n==((Ij[n2][n])-1))&&(n4==((Ji[n][n4])-1))){						//J在I的覆蓋範圍內，J'在I的覆蓋範圍內
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


				if ((n==((Ij[n2][n])-1))){																		//在i的覆蓋範圍內
						

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
	    Coes.add(Pn/(4.11*pow(10.0,(-14.0))) );																	//Pn = 4.11*pow(10.0,(-14.0)	//定義在上方define
		Vars.add(vars[  M*N + n2 ]);																			//SINR*u

	//   cout << "(Pn/(2*pow(10.0,(-14))) )=" << (Pn/(2*pow(10.0,(-14))) ) <<endl;
	//Vars.add(vars[  M*N + n*N + n2 ]);   //sinr
		cout << "SINR_D_%d" << n2+1 <<endl;
		sprintf(cst_name, "SINR_D_%d",n2+1);
		IloConstraint cst;
		cst = IloScalProd(Coes,Vars) == 0;																		//限制式
		cst.setName(cst_name);
		model.add(cst);
		Vars.clear();
		Coes.clear();

		



	node2 = node2-> GetNext();
}  

//[8] Flow 根據Shannon capacity ，我們可以計算BS i和TP j 之通道容量。

   node2  = this-> GetNode_T();												//用node的位址找到該TP
		while (node2 != NULL){  
			n2 = node2->GetValue()->GetNo();								//j
 
				Coes.add(-W / (log(2.0)));									//w/ln(2)
				Vars.add(vars[ M*N + N + n2 ]);								//H*u的記憶體位置

				Coes.add(Rate);												//R
				Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);						//Zj*u

				char cst_name[60];
				sprintf(cst_name, "Capacity_%d ", n2+1  );             
				IloConstraint cst;                                      //創建限制式cst
				cst = IloScalProd(Coes,Vars) <= 0;                      //生成一條限制式
				cst.setName(cst_name);                                  
				model.add(cst);                                         //將限制式加入模型中
			    Coes.clear();                                           //清空Coes陣列
                Vars.clear();                                           //清空Vars陣列
				  
	    	node2 = node2-> GetNext();
	      }   
	   


//-------------------------------------------------------Biliner Tearm : V--------------------------------------------

//[9.1] λv*u
 
node  = this-> GetNode_B();                                             //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                                      //i
	 
    node4  = this-> GetNode_T();                                        //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                //j'
 		 
				for(int j=0 ; j<par ; j++)
				{
					Coes.add(1);
				    Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //新增λv*u的變數
                } 

				Coes.add(-1);
				Vars.add(vars[ 2*M*N + 3*N + M ]);                                  //新增u的變數
				
				char cst_name[60];
				sprintf(cst_name, "Sumλv_%d_%d ",n+1, n4+1 );          //把Sumλv_%d.%d寫入cst_name
				IloConstraint cst;                                      //創建限制式cst
				cst = IloScalProd(Coes,Vars) ==0;                       //生成一條限制式
				cst.setName(cst_name);                                  //設立限制式名稱
				model.add(cst);                                         //將限制式加入模型中
			    Coes.clear();                                           //清空Coes陣列
                Vars.clear();                                           //清空Vars陣列
				  
		  
	    	node4 = node4-> GetNext();
	       }   
	   
    	    node = node -> GetNext();                                       
	       }
           //總共M*N條限制式
 
//[9.2.1]
			 
node2  = this-> GetNode_T();                                                       //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                               //j
	 
 	   node = this-> GetNode_B();                                                  //用node的位址找到該BS
		   while(node != NULL){
		   n = node->GetValue()->GetNo();                                          //i
				
		      node4 = this-> GetNode_T();                                          //用node的位址找到該TP
		            while(node4 != NULL){
			        n4 = node4->GetValue()->GetNo();                               //j'
			
			if (n2!=n4)                                                            // j!=j' 
			{                                                 
				for(int j=0 ; j<par ; j++)
			    {				 									
					Coes.add(-1);
	         		Vars.add (bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);					//新增△V*u par的變數 
			    }  
		
				    Coes.add(1);
					Vars.add (vars[ M*N + n2 ]);                                                        //新增SINR*u的變數
	
					sprintf(cst_name, "Sum_△V*u_%d_%d_%d",n2+1,n+1,n4+1);								//把Sum_△V*u_%d_%d_%d寫入cst_name
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
			 
		   	//總共N*M*N條限制式
		    
 
//[9.2.2] SINR_△V
			 
node2  = this-> GetNode_T();                                                    //用node的位址找到該TP
	 while(node2 != NULL){
	 n2 = node2->GetValue()->GetNo();                                           //j
	 
		node = this-> GetNode_B();                                              //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                                      //i
				 

		node4 = this-> GetNode_T();                                             //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();                                    //j'
		 
				    
			 if (n2!=n4)                                                        // j!=j' 
			 {                                               
				 for(int j=0 ; j<par ; j++)
				 {
 
					Coes.add(1);
	         		Vars.add (bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);					//新增△V*u par的變數
		             
	          
					Coes.add(-(Var_UB [ M*N + n2 ] - Var_LB[ M*N + n2 ]));                              // (SINR_UB - SINR_LB)*u
					Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);						//新增λv*u par的變數
 
                    char cst_name[60];
			        sprintf(cst_name,  "△S_V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);							//把△S_V_%d_%d_%d_%d寫入cst_name
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
			 	             
            //總共N*M*N條限制式
 

//RLTV     
//[9.3.1] RLTV1		
 		 
 node2  = this-> GetNode_T();                                                 //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                          //j
	 	 
		node = this-> GetNode_B();                                            //用node的位址找到該BS
		   while(node != NULL){
		   n = node->GetValue()->GetNo();                                     //i
				 
		       node4 = this-> GetNode_T();                                    //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);							//新增△V*u par的變數
			}
	
					Coes.add(-1*Var_LB[ M*N + n2 ]);															//SINR_LB
					Vars.add(vars[ n*N + n4 ]);																	//新增Power*u的變數

					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);													//新增V*u的變數
                    
                    sprintf(cst_name, "V:RLT1_%d_%d_%d",n2+1,n+1,n4+1);											//把V:RLT1_%d_%d_%d寫入cst_name
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
			 //總共N*M*N條限制式

  
//[9.3.2] RLTV2   

node2  = this-> GetNode_T();                                        //用node的位址找到該TP
	 while(node2 != NULL){
     n2 = node2->GetValue()->GetNo();                               //j
	 		 
		node = this-> GetNode_B();                                  //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                          //i
				 
		       node4 = this-> GetNode_T();                          //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);					//新增△V*u par的變數
		    }
	
                  	Coes.add( -1*Var_LB[ M*N + n2 ]);                                                   //SINR_LB
					Vars.add(vars[ n*N + n4 ]);                                                         //新增Power*u的變數
            
					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);                                            //新增V*u的變數
				 
					char cst_name[60];
					sprintf(cst_name, "V:RLT2_%d_%d_%d", n2+1,n+1,n4+1);                                //把V:RLT2_%d_%d_%d寫入cst_name
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
	   //總共N*M*N條限制式
 
//[9.3.3] RLTV3 
	 
node2  = this-> GetNode_T();                               //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                       //j
	 	 
		node = this-> GetNode_B();                         //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                 //i
				 
		       node4 = this-> GetNode_T();                 //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);						//新增△V*u par的變數

					Coes.add ( Interval*(Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ]) );                        //Interval*(SINR_UB-SINR_LB)
					Vars.add (vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);                          //新增λv*u par的變數 
		 
         	  }	
		  					
					Coes.add(-1*Var_UB[ M*N + n2 ]);                                                        //SINR_UB
					Vars.add(vars[ n*N + n4 ]);                                                             //新增Power*u的變數

					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);												//新增V*u的變數
				    
					sprintf(cst_name, "V:RLT3_%d_%d_%d",n2+1,n+1,n4+1 );                                    //把V:RLT3_%d_%d_%d寫入cst_name
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
	   //總共N*M*N條限制式
	
  
//[9.3.4] RLTV4
 		 
node2  = this-> GetNode_T();                               //用node的位址找到該TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                       //j
	 
 	   node = this-> GetNode_B();                          //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();                 //i
				 
		      node4 = this-> GetNode_T();                  //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);							//新增△V*u par的變數

					Coes.add ( Interval*(Var_UB[ M*N + n2 ]- Var_LB[ M*N + n2 ]));								//Interval*(SINR_UB-SINR_LB)
					Vars.add ( vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);								//新增λv*u par的變數
		 
         	  }	
		  

					Coes.add(-1*Var_UB[ M*N + n2 ]);                                                            //SINR_UB
					Vars.add(vars[ n*N + n4 ]);                                                                 //新增Power(i,j')*u的變數

					Coes.add(1);
	         		Vars.add (bilvars[ n2*M*N + n*N + n4 ]);                                                    //新增V*u的變數
				    
					sprintf(cst_name, "V:RLT4_%d_%d_%d",n2+1,n+1, n4+1 );										//把V:RLT4_%d_%d_%d寫入cst_name
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
	   //總共N*M*N條限制式 
	
    
//-------------------------------------------------------Bilinear Tearm :K---------------------------------------//

//[10.1] λk*u
 
node3  = this-> GetNode_B();                                             //用node的位址找到該BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();                                     //i'
	 
    node4  = this-> GetNode_T();                                         //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();                                 //j'
 		 
				for(int j=0 ; j<par ; j++)
				{
					Coes.add(1);
				    Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);   //新增λk*u的變數
                } 

				Coes.add(-1);
				Vars.add(vars[ 2*M*N + 3*N + M ]);                                             //新增u的變數
				
				char cst_name[60];
				sprintf(cst_name, "Sumλk_%d_%d ",n3+1, n4+1 );         //把Sumλk_%d.%d寫入cst_name
				IloConstraint cst;                                      //創建限制式cst
				cst = IloScalProd(Coes,Vars) ==0.0;                     //生成一條限制式
				cst.setName(cst_name);                                  //設立限制式名稱Sumλ_%d_%d
				model.add(cst);                                         //將限制式加入模型中
			    Coes.clear();                                           //清空Coes陣列
                Vars.clear();                                           //清空Vars陣列
				  
		  
	    	node4 = node4 -> GetNext();
	       }   
	   
    	    node3 = node3 -> GetNext();                                       
	       }
           //總共M*N條限制式


//[10.2.1]
			  
node2  = this-> GetNode_T();                                     //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                             //j
	 	 
		node3 = this-> GetNode_B();                              //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                     //i'
				 
		       node4 = this-> GetNode_T();                       //用node的位址找到該TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();              //j'
			
			if (n2!=n4)        // j!=j' 
			{                       
				for(int j=0 ; j<par ; j++)
				{					
					
					Coes.add(-1);    
	         		Vars.add (bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);       //新增△K*u的變數	  
		        
				}  
                
				    Coes.add(1);
					Vars.add (vars[ M*N + n2 ]);									//新增SINR*u的變數 
                    	
					sprintf(cst_name, "Sum_△K*u_%d_%d_%d",n2+1,n3+1,n4+1);			//把Sum_△K*u_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式
		    

//[10.2.2] SINR_△K   
 		 
node2  = this-> GetNode_T();                             //用node的位址找到該TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                     //j
	 	
		node3 = this-> GetNode_B();                      //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();             //i'
				 
		       node4 = this-> GetNode_T();               //用node的位址找到該TP
		           while(node4 != NULL){
			       n4 = node4->GetValue()->GetNo();      //j'
		 

		if (n2!=n4)      // j!=j'
		{                         
			for(int j=0 ; j<par ; j++)
			{
										
				    Coes.add(1);                           
	         		Vars.add (bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);			//新增△K*u的變數
					                 
				    Coes.add(-(Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ]));
				    Vars.add( vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);				//新增λk*u的變數

					sprintf(cst_name,  "△S_K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                             //把△S_K_%d_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式
  
//RLTK
//[10.3.1] RLTK1  		
 	 
node2  = this-> GetNode_T();                                              //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                      //j
	 
 		 
	   node3 = this-> GetNode_B();                                        //用node的位址找到該BS
		   while(node3 != NULL){
		   n3 = node3->GetValue()->GetNo();                               //i'
				 

		     node4 = this-> GetNode_T();                                  //用node的位址找到該TP
		         while(node4 != NULL){
			     n4 = node4->GetValue()->GetNo();                         //j'
	 

		if (n2!=n4)    // j!=j'
		{                   

			for(int j=0 ; j<par ; j++)
			{
				    double a1 = (j+1)/nsum;
				    double a2 = j/nsum;
				    double a = (pow(a1,ga)-pow(a2,ga))*(Var_UB[ n3*N + n4 ] - Var_LB[ n3*N + n4 ]);		//方法1  *


					Coes.add( -1*(Var_LB[ n3*N + n2 ] + a*j) );                                         //Power_LB + a*(np-1)
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);		//新增△K*u的變數

			}	
				    
			        Coes.add(-1*Var_LB[ M*N + n2 ]);													//SINR_LB
					Vars.add(vars[ n3*N + n4 ]);														//新增Power*u的變數

					Coes.add(1);
  					Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);									//新增K*u的變數
			 
                    sprintf(cst_name, "K:RLT1_%d_%d_%d", n2+1,n3+1,n4+1);                               //把K:RLT1_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式 
		   	  
   
//[10.3.2] RLTK2
			 
node2  = this-> GetNode_T();                                      //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                              //j
	 
		node3 = this-> GetNode_B();                               //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                      //i'
				 
               node4 = this-> GetNode_T();                        //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j  ]);     //新增△K*u的變數
		 
		    }	
	
			        Coes.add(-1*Var_LB[ M*N + n2]);														//SINR_LB
					Vars.add(vars[ n3*N + n4 ]);														//新增Power*u的變數

					Coes.add(1);
  					Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);                                   //新增K*u的變數		

                    sprintf(cst_name, "K:RLT2_%d_%d_%d",n2+1,n3+1,n4+1);								//把K:RLT2_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式 
	

//[10.3.3] RLTK3 
		 
node2  = this-> GetNode_T();                                          //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                                  //j
	 
 		 
		node3 = this-> GetNode_B();                                   //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                          //i'
				 

		       node4 = this-> GetNode_T();                            //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);		//新增△K*u的變數
			 

					Coes.add ( Interval*(Var_UB[ M*N + n2]- Var_LB[ M*N + n2]));						//Interval*(SINR_UB-SINR_LB)
					Vars.add ( vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);               //新增λk*u的變數
		 
              }	
						
				    Coes.add(-1*Var_UB[ M*N + n2 ]);													//SINR_UB
					Vars.add(vars[ n3*N + n4 ]);														//新增Power*u的變數

					Coes.add(1);
	         		Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);									//新增K*u的變數
				    
					sprintf(cst_name, "K:RLT3_%d_%d_%d", n2+1,n3+1,n4+1 );								//把K:RLT3_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式 
	
  
//[10.3.4] RLTK4
			 
node2  = this-> GetNode_T();                                   //用node的位址找到該TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                           //j
	 		 
		node3 = this-> GetNode_B();                            //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();                   //i'
				 
		       node4 = this-> GetNode_T();                     //用node的位址找到該TP
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
					Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);		//新增△K*u的變數
		 

					Coes.add ( Interval*(Var_UB[ M*N + n2 ]- Var_LB[ M*N + n2 ]));                      //Interval*(SINR_UB-SINR_LB)
					Vars.add ( vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);          //新增λk*u的變數
		 
         	  }	
		  				
				    Coes.add(-1*Var_UB[ M*N + n2]);														//SINR_UB
					Vars.add(vars[ n3*N + n4 ]);														//新增Power*u的變數

					Coes.add(1);
	         		Vars.add (bilvars[ N*M*N + n2*M*N + n3*N + n4 ]);									//新增K*u的變數
				 
					sprintf(cst_name, "K:RLT4_%d_%d_%d",n2+1,n3+1, n4+1 );								//把K:RLT4_%d_%d_%d寫入cst_name
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
           //總共N*M*N條限制式
 

//[12.2 - 12.5] Capacity (包含移除對數項目，要記得看求解方式)
 
    float HUB, HLB, slmMinus,R , D, Q , Padd1;	
 

//[12.2]
			 
node2  = this-> GetNode_T();										//用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();								//j
	 
 		 
				
					 Padd1 = 1 + Var_LB[ M*N + n2 ];                //1+SINR_LB
					 R = log(Padd1)-1;
					 Q = R*Padd1 + 1;
					
					//cout<< "Pmax =" <<Pmax<<endl;

					Coes.add(Padd1);                                //1+SINR_LB
					Vars.add(vars[ M*N + N + n2 ]);					//新增H*u的變數
					
					Coes.add(-1);
					Vars.add(vars[ M*N + n2 ]);						//新增SINR*u的變數
					
					Coes.add(-Q);
					Vars.add(vars[ 2*M*N + 3*N + M ]);              //新增u的變數
										
                    
					char cst_name[60];
					sprintf(cst_name, "Capacity1_D_%d",n2+1);		//把Capacity1_D_%d寫入cst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0 ;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();

			 node2 = node2 -> GetNext();	
		    }
			 
            
	        //總共N條限制式

//[12.3]
		 
node2  = this-> GetNode_T();             //用node的位址找到該TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();     //j
	 
				    HUB = 1 + Var_UB[ M*N + n2 ];							//HUB=1+SINR_UB
				    HLB = 1 + Var_LB[ M*N + n2 ];							//HLB=1+SINR_LB
				    slmMinus = Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ];		//SINR_UB-SINR_LB
			
 
				    Padd1 = 1 + beta;										//beta為目標式				
					R = log(Padd1)-1;
					Q = R*Padd1 + 1;
 
					Coes.add(Padd1);
					Vars.add(vars[ M*N + N + n2 ]);							//新增H*u的變數
					
					Coes.add(-1);
					Vars.add(vars[ M*N + n2 ]);								//新增SINR*u的變數
					
					Coes.add(-Q);
					Vars.add(vars[ 2*M*N + 3*N + M ]);						//新增u的變數
				
                    char cst_name[60];
					sprintf(cst_name, "Capacity2_%d", n2+1);				//把Capacity2_%d寫入cst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();
				
			 node2 = node2 -> GetNext();	
		    }
	
	        //總共N條限制式
					
//[12.4]
			 
node2  = this-> GetNode_T();												//用node的位址找到該TP  
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();										//j

					Padd1 = 1 + Var_UB[ M*N + n2 ];							//Padd1=1+SINR_UB
					R = log(Padd1)-1;
					Q = R*Padd1 + 1;
										
					Coes.add(Padd1);
					Vars.add(vars[ M*N + N + n2 ]);							//新增H*u的變數
					
					Coes.add(-1);
					Vars.add(vars[ M*N + n2 ]);								//新增SINR*u的變數

					Coes.add(-Q);
					Vars.add(vars[ 2*M*N + 3*N + M ]);						//新增u的變數
 
                    char cst_name[60];
					sprintf(cst_name, "Capacity3_%d", n2+1);				//把Capacity3_%d寫入cst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) <= 0;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();
					
			 node2 = node2 -> GetNext();	
		    } 
			 
	        //總共N條限制式
					
//[12.5] 
			 
node2  = this-> GetNode_T();                 //用node的位址找到該TP
	while(node2 != NULL){  
	n2 = node2->GetValue()->GetNo();         //j
			      
					slmMinus = Var_UB[ M*N + n2 ] - Var_LB[ M*N + n2 ];         //SINR_UB-SINR_LB
		     		HUB = log(1 + Var_UB[ M*N + n2 ] );                         //ln(1+SINR_UB)
					HLB = log(1 + Var_LB[ M*N + n2 ] );                         //ln(1+SINR_LB)
					R =  Var_UB[ M*N + n2 ]*HLB;                                //SINR_UB*HLB
					D =  Var_LB[ M*N + n2 ]*HUB;                                //SINR_LB*HUB
            

					Coes.add(slmMinus);
					Vars.add(vars[ M*N + N + n2 ]);								//新增H*u的變數

					Coes.add(HLB-HUB);
					Vars.add(vars[ M*N + n2 ]);									//新增SINR*u的變數

					Coes.add(D-R);
					Vars.add(vars[ 2*M*N + 3*N + M ]);							//新增u的變數

                    char cst_name[60];
		            sprintf(cst_name, "Capacity4_%d", n2+1);					//把Capacity4_%d寫入cst_name
					IloConstraint cst;
					cst = IloScalProd(Coes,Vars) >= 0 ;
					cst.setName(cst_name);
					model.add(cst);
					Coes.clear();
					Vars.clear();
				
		     node2 = node2 -> GetNext();	
		    }
			
            //總共N條限制式

//乘u限制式
//[13.1] y*u

 node  = this-> GetNode_B();                   //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();             //i
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + n ]);       //新增y*u的變數

			Coes.add(-1);                      
			Vars.add(vars[ 2*M*N + 3*N + M ]);       //新增u的變數

			sprintf(cst_name, "y_u_%d",n+1);         //把y_u_%d寫入cst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node  = node  -> GetNext();
	    }
	    //總共M條限制式 

//[13.2] x*u

 node  = this-> GetNode_B();               //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();         //i
     	 	    
		node2  = this-> GetNode_T();       //用node的位址找到該TP
		while (node2 != NULL){  
		n2 = node2->GetValue()->GetNo();   //j 
		    	

			Coes.add(1);
			Vars.add(vars[ M*N + 2*N + n*N + n2 ]);        //新增x*u的變數

			Coes.add(-1); 
			Vars.add(vars[ 2*M*N + 3*N + M ]);             //新增u的變數

			sprintf(cst_name, "x_u_%d_%d",n+1,n2+1);       //把x_u_%d_%d寫入cst_name
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
	       //總共M*N條限制式

//[13.3] Zj*u

 node2  = this-> GetNode_T();                   //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();            //j
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);   //新增Zj*u的變數

			Coes.add(-1);                      
			Vars.add(vars[ 2*M*N + 3*N + M ]);        //新增u的變數
            
			sprintf(cst_name, "Zj_u_%d",n2+1);        //把Zj_u_%d寫入cst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node2  = node2  -> GetNext();
	    }
	    //總共N條限制式
 
//[13.4.1] λv*u

 node  = this-> GetNode_B();                      //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                //i
     	 	    
		node4  = this-> GetNode_T();              //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();          //j' 
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //新增λv*u的變數

			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);                              //新增u的變數
		
			sprintf(cst_name, "λv_u_%d_%d_%d",n+1,n4+1,j+1);               //把λv_u_%d_%d_%d寫入cst_name
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
	      //總共M*N條限制式

//[13.4.2] λk*u

 node3 = this-> GetNode_B();                      //用node的位址找到該BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();              //i'
     	 	    
		node4  = this-> GetNode_T();              //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();          //j' 
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);   //新增λk*u的變數

			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);                                         //新增u的變數
		
			sprintf(cst_name, "λk_u_%d_%d_%d",n3+1,n4+1,j+1);                         //把λk_u_%d_%d_%d寫入cst_name
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
	      //總共M*N條限制式


//[13.5] △V*u
/*		 
node2  = this-> GetNode_T();                     //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();             //j
	 		 
	    node = this-> GetNode_B();               //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();       //i
				 
		node4 = this-> GetNode_T();              //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();     //j'
        
        
			for(int j=0 ; j<par ; j++)
			{
				Coes.add(1);
			    Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);   //新增△V*u的變數

			    Coes.add(-1);
			    Vars.add(vars[ 2*M*N + 3*N + M ]);                                                //新增u的變數
		
			    sprintf(cst_name, "△V_u_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1 );                        //把△V_u_%d_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式
 
//[13.6] △K*u
			 
node2  = this-> GetNode_T();                   //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();           //j
	 
		node3 = this-> GetNode_B();            //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();   //i'
				 
		node4 = this-> GetNode_T();            //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();   //j'
		 
        
			for(int j=0 ; j<par ; j++)
			{
				Coes.add(1);
			    Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);    //新增△K*u的變數

			    Coes.add(-1);
		      	Vars.add(vars[ 2*M*N + 3*N + M ]);                                                  //新增u的變數
		
			    sprintf(cst_name, "△K_u_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1 );                         //把△K_u_%d_%d_%d_%d寫入cst_name
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
			//總共N*M*N*par條限制式

*/
 
//[14] Big M-------------------------------------------------------------------------
 
//[14.1] y*u

 node  = this-> GetNode_B();                      //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                //i
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + n ]);          //新增y*u的變數

			Coes.add(-Big);                             //Big M
			Vars.add(binvars[ n ]);                     //新增y的變數
            
			sprintf(cst_name, "Big M_y_%d",n+1);        //把Big M_y_%d寫入cst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

	    node  = node  -> GetNext();
	   }
	   //總共M條限制式
 
//[14.2] x*u

 node  = this-> GetNode_B();                       //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                 //i
     	 	    
		node2  = this-> GetNode_T();               //用node的位址找到該TP
		while (node2 != NULL){  
		n2 = node2->GetValue()->GetNo();           //j
		    	

			Coes.add(1);
			Vars.add(vars[ M*N + 2*N + n*N + n2 ]);    //新增x*u的變數

			Coes.add(-Big);                            //Big M
			Vars.add(binvars[ M + n*N + n2 ]);         //新增x的變數

			sprintf(cst_name, "Big M_x_%d_%d",n+1,n2+1);   //把Big M_x_%d_%d寫入cst_name
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
	       //總共M*N條限制式

//[14.3] Zj*u

 node2 = this-> GetNode_T();                      //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();              //j
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);       //新增Zj*u的變數

			Coes.add(-Big);                         //Big M
			Vars.add(binvars[ M + M*N + n2 ]);      //新增Zj的變數

			sprintf(cst_name, "Big M_Zj_%d",n2+1);  //把Big M_Zj_%d寫入cst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= 0;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

	    node2 = node2 -> GetNext();
	   }
	   //總共N條限制式
 
//[14.4.1] λv*u

 node  = this-> GetNode_B();                        //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                  //i
     	 	    
		node4  = this-> GetNode_T();                //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();            //j'
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //新增λv*u的變數

			Coes.add(-Big);                                           //Big M
			Vars.add(binvars[ M + M*N + N + n*N*par + n4*par + j ]);  //新增λv的變數
		
			sprintf(cst_name, "Big M_λv_%d_%d_%d",n+1,n4+1,j+1);     //把Big M_λv_%d_%d_%d寫入cst_name
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
	      //總共M*N條限制式
  
//[14.4.2] λk*u

 node3 = this-> GetNode_B();                        //用node的位址找到該BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();                //i'
     	 	    
		node4  = this-> GetNode_T();                //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();            //j'
		    	
		for(int j=0 ; j<par ; j++)
		{
			Coes.add(1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);   //新增λk*u的變數

			Coes.add(-Big);                                                      //Big M
			Vars.add(binvars[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ]);  //新增λk的變數
		
			sprintf(cst_name, "Big M_λk_%d_%d_%d",n3+1,n4+1,j+1);    //把Big M_λk_%d_%d_%d寫入cst_name
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
	      //總共M*N條限制式
  

//[14.5] △V*u  
/*		 
node2  = this-> GetNode_T();                       //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();               //j
 		 
		node = this-> GetNode_B();                 //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();         //i
				 
		node4 = this-> GetNode_T();                //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();       //j'
		   
		if(n2!=n4) // j!=j'
		{                            
			for(int j=0 ; j<par ; j++)
		    {
				Coes.add(1);
			    Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);   //新增△V*u的變數

			    Coes.add(-Big);                                                                   //Big M
			    Vars.add(bilvars[ 2*N*M*N + 2*N*M*N*par + n2*M*N*par + n*N*par + n4*par + j ]);                 //新增

				sprintf(cst_name, "Big M_△V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);                     //把Big M_△V_%d_%d_%d_%d寫入cst_name
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
	        //總共N*M*N條限制式
 

//[14.6] △K*u 
	 
	node2  = this-> GetNode_T();              //用node的位址找到該TP
	    while(node2 != NULL){
		n2 = node2->GetValue()->GetNo();      //j
 		 
		node3 = this-> GetNode_B();           //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();  //i'
				 
		node4 = this-> GetNode_T();           //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();  //j'

			if(n2!=n4)    // j!=j'
			{
				for(int j=0 ; j<par ; j++)
		        {
					Coes.add(1);
			        Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);      //新增△K*u的變數

			        Coes.add(-Big);                                                                       //Big M
			        Vars.add(bilvars[ 2*N*M*N + 3*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);        //新增
 
					sprintf(cst_name, "Big M_△K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                        //把Big M_△K_%d_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式
		   

 */
//[15] Big M 
 
//[15.1] y*u

 node  = this-> GetNode_B();                       //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();                 //i
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);           //新增u的變數

			Coes.add(Big);                         //Big M
			Vars.add(binvars[ n ]);                //新增y的變數
			
			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 2*N + n ]);     //新增y*u的變數
            
			sprintf(cst_name, "Mix M_y_%d",n+1);   //把Mix M_y_%d寫入cst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node  = node  -> GetNext();
	   }
       //總共M條限制式


//[15.2] x*u

 node  = this-> GetNode_B();                     //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();               //i
     	 	    
		node2  = this-> GetNode_T();             //用node的位址找到該TP
		while (node2 != NULL){  
		n2 = node2->GetValue()->GetNo();         //j
		    	

		    Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);          //新增u的變數

			Coes.add(Big);                        //Big M
			Vars.add(binvars[ M + n*N + n2 ]);    //新增x的變數
			
			Coes.add(-1);
     		Vars.add(vars[ M*N + 2*N + n*N + n2 ]);   //新增x*u的變數
 
			sprintf(cst_name, "Mix_M_x_%d_%d",n+1,n2+1);  //把Mix_M_x_%d_%d寫入cst_name
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
	       //總共M*N條限制式


//[15.3] Zj*u

 node2  = this-> GetNode_T();                        //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();                 //j
        	
            Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);             //新增u的變數

			Coes.add(Big);                           //Big M
			Vars.add(binvars[ M + M*N + n2 ]);       //新增的Zj變數
			
			Coes.add(-1);
			Vars.add(vars[ 2*M*N + 2*N + M + n2 ]);        //新增Zj*u的變數
            
			sprintf(cst_name, "Mix M_Zj_%d",n2+1);   //把Mix M_Zj_%d寫入cst_name
 			IloConstraint cst;
			cst = IloScalProd(Coes,Vars) <= Big;
			cst.setName(cst_name);
			model.add(cst); 
		   	Coes.clear();
			Vars.clear();	

		 node2  = node2 -> GetNext();
	   }
       //總共N條限制式

 
//[15.4.1] λv*u

 node  = this-> GetNode_B();                  //用node的位址找到該BS
	while(node != NULL){
	n = node->GetValue()->GetNo();            //i
     	 	    
		node4  = this-> GetNode_T();          //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();      //j'
		    	
		for(int j=0 ; j<par ; j++)
		{

			Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);								//新增u的變數

			Coes.add(-1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + n*N*par + n4*par + j ]);   //新增λv*u的變數
 
			Coes.add(Big);													//Big M
			Vars.add(binvars[ M + M*N + N + n*N*par + n4*par + j ]);		//新增λv的變數
		
			sprintf(cst_name, "Mix_M_λv_%d_%d_%d",n+1,n4+1,j+1);			//把Mix_M_λv_%d_%d_%d寫入cst_name
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
           //總共M*N條限制式


//[15.4.2] λk*u

 node3  = this-> GetNode_B();                  //用node的位址找到該BS
	while(node3 != NULL){
	n3 = node3->GetValue()->GetNo();           //i'
     	 	    
		node4  = this-> GetNode_T();           //用node的位址找到該TP
		while (node4 != NULL){  
		n4 = node4->GetValue()->GetNo();       //j'
		    	
		for(int j=0 ; j<par ; j++)
		{

			Coes.add(1);
			Vars.add(vars[ 2*M*N + 3*N + M ]);      //新增u的變數

			Coes.add(-1);
			Vars.add(vars[ 1 + 2*M*N + 3*N + M + M*N*par + n3*N*par + n4*par + j ]);    //新增λk*u的變數
 
			Coes.add(Big);                                                        //Big M
			Vars.add(binvars[ M + M*N + N + M*N*par + n3*N*par + n4*par + j ]);   //新增λk的變數
		
			sprintf(cst_name, "Mix_M_λk_%d_%d_%d",n3+1,n4+1,j+1);                //把Mix_M_λk_%d_%d_%d寫入cst_name
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
           //總共M*N條限制式


//[15.5] △V*u
/*	 
node2  = this-> GetNode_T();                    //用node的位址找到該TP
	while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();            //j
	 
		node = this-> GetNode_B();              //用node的位址找到該BS
			while(node != NULL){
			n = node->GetValue()->GetNo();      //i
				 
		node4 = this-> GetNode_T();             //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();    //j'
		    	
			if(n2!=n4)  // j!=j'
			{    	                
				for(int j=0 ; j<par ; j++)
				{
					Coes.add(1);
			        Vars.add(vars[ 2*M*N + 3*N + M ]);        //新增u的變數
			
			        Coes.add(-1);
			        Vars.add(bilvars[ 2*N*M*N + n2*M*N*par + n*N*par + n4*par + j ]);    //新增△V*u的變數

			        Coes.add(Big);                                                                     //Big M
			        Vars.add(bilvars[ 2*N*M*N + 2*N*M*N*par + n2*M*N*par + n*N*par + n4*par + j ]);             //新增

			        sprintf(cst_name, "Mix_M_△V_%d_%d_%d_%d",n2+1,n+1,n4+1,j+1);                      //把Mix_M_△V_%d_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式

 
//[15.6] △K*u
		 
node2  = this-> GetNode_T();                    //用node的位址找到該TP
    while(node2 != NULL){
	n2 = node2->GetValue()->GetNo();            //j
 		 
		node3 = this-> GetNode_B();             //用node的位址找到該BS
			while(node3 != NULL){
			n3 = node3->GetValue()->GetNo();    //i'
				 
		node4 = this-> GetNode_T();             //用node的位址找到該TP
		    while(node4 != NULL){
			n4 = node4->GetValue()->GetNo();    //j'
		    	
		if(n2!=n4)    // j!=j' 
		{
			for(int j=0 ; j<par ; j++)
		    {
				Coes.add(1);
			    Vars.add(vars[ 2*M*N + 3*N + M ]);        //新增u的變數
			
			    Coes.add(-1);
			    Vars.add(bilvars[ 2*N*M*N + N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);   //新增△K*u的變數

			    Coes.add(Big);                                                                     //Big M
			    Vars.add(bilvars[ 2*N*M*N + 3*N*M*N*par + n2*M*N*par + n3*N*par + n4*par + j ]);     //新增

			    sprintf(cst_name, "Mix_M_△K_%d_%d_%d_%d",n2+1,n3+1,n4+1,j+1);                     //把Mix_M_△K_%d_%d_%d_%d寫入cst_name
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
			//總共N*M*N條限制式
*/		   
 
//--------------------------------------------------------------------------------
 


IloCplex cplex(model);
//cplex.setParam(cplex.EpAGap,0.05);
cplex.solve();
//cplex.setParam(cplex.TiLim,10000); //限制跑的時間
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
 

//顯示值------------------------------------------------------------------------------------------------------
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

//Bin λv
	    	
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

//Bin λk      
		
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

//λv*u

		    fprintf(fp,"[λv_u]\n");	
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

//λk*u

            fprintf(fp,"[λk_u]\n");	
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


fprintf(fp,"[△V]\n");	
	 		 
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
fprintf(fp,"[△K]\n");		 
			 
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

