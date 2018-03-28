#pragma once 
#include <iostream>
#include <ilcplex/ilocplex.h>

template <class TYPE>
class LIST{
	private:
		LIST *Prev;			//�e�@�Ӥ���
		LIST *Next;			//��@�Ӥ���
		TYPE Value;			//��
	public:
		LIST(){				//�غc�l�G	
			Prev = NULL;
			Next = NULL;
		}

		LIST(TYPE value){	//�غc�l�G��
			Prev = NULL;
			Next = NULL;
			Value = value;
		}

		LIST(LIST<TYPE> *prev, TYPE value, LIST<TYPE> *next){	//�غc�l�G�e����,��,�ᤸ��
			Prev = prev;
			Next = next;
			Value = value;
		}

		~LIST(){			//�Ѻc�l
			if(Prev != NULL)
				Prev->SetNext(Next);
			if(Next != NULL)
				Next->SetPrev(Prev);
		}

		TYPE GetValue(void){return Value;}	//���o�G��
			
		TYPE *GetValueRef(void){return &Value;}	//���o�G��
		
		LIST *GetPrev(void){return Prev;}	//���o�G�e�@�Ӥ���
		
		LIST *GetNext(void){return Next;}	//���o�G��@�Ӥ���

		void SetValue(TYPE value){Value = value;}	//�]�w�G��
		
		void SetPrev(LIST *prev){Prev = prev;}	//�]�w�G�e�@�Ӥ���
		
		void SetNext(LIST *next){Next = next;}	//�]�w�G��@�Ӥ���
		
		LIST *Head(void){		//���o�G��������
			LIST *current=this;
			while(current->GetPrev() != NULL)
				current = current->GetPrev();
			return current;
		}

		LIST *Tail(void){		//���o�G��������
			LIST *current=this;
			while(current->GetNext() != NULL)
				current = current->GetNext();		
			return current;
		}

		LIST *Find(TYPE value){	//�M��G�ȩҦb���̫e�ݤ���
			LIST<TYPE> *current=Head();
			do{
				if(current->GetValue() == value)
					return current;
				else;
				current = current->GetNext();
			}while(current != NULL);
			return NULL;
		}

		int GetLength(void){		//���o�GLIST����(�����ƶq)
			int i=0;
			LIST<TYPE> *current;
			current = Head();
			while(current!=NULL){
				i++;
				current = current->GetNext();
			}
			return i;
		}

		void Remove(){			//�����G�Ӷ�����
			if(Prev != NULL)
				Prev->SetNext(Next);
			if(Next != NULL)
				Next->SetPrev(Prev);
			//cout<<"asdfasdfsdaf"<<endl;
			delete this;
		}

		LIST& operator+(LIST& RHS){	//�B��l�����G + �s����LIST
			RHS.Head()->SetPrev(this->Tail());
			Tail()->SetNext(RHS.Head());
		}
};


