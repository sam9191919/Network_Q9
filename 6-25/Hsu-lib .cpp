#pragma once 
#include "lib.h"


/*
LIST& LIST<TYPE> :: operator+(LIST& RHS){	//�B��l�����G + �s����LIST
	RHS.Head()->SetPrev(this->Tail());
	Tail()->SetNext(RHS.Head());
}
*/