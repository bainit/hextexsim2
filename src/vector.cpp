#include <iostream>
#include <utility>
#include <memory>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_util.h"

using namespace std;

Vector::Vector(Orientation o) : orientation{o}  {

	vector = create_vector();

}

Vector::Vector(initializer_list<double> arr, Orientation o) : orientation{o} {

	vector = create_vector();
	
	initializer_list<double>::iterator itr = arr.begin();

	if(ishorizontal()) {	

		int i=0;		
		try{
			for(int j=0;j<3;j++) {
				gsl_matrix_set(vector,i,j,*itr);
				itr++;
			}
		}catch(runtime_error e) {
			cout<<e.what()<<endl;
		}

		if(itr!=arr.end()) {
			gsl_matrix_free(vector);
			throw runtime_error("Vector::Vector(initializer_list): to many elements in an array");

		}

	}else{
		int i=0;
		try{
			for(int j=0;j<3;j++) {
				gsl_matrix_set(vector,j,i,*itr);
				itr++;
			}
		}catch(runtime_error e) {
			cout<<e.what()<<endl;
		}

		if(itr!=arr.end()) {
			gsl_matrix_free(vector);
			throw runtime_error("Vector::Vector(initializer_list): to many elements in an array");

		}

	}

}

Vector::Vector(gsl_matrix *vect, Orientation o) : orientation{o} {

	vector = create_vector();
		
	if(vect->size1==vector->size1 and vect->size2==vector->size2) {

		gsl_matrix_memcpy(vector,vect);
	
	}else {
		gsl_matrix_free(vector);
		throw runtime_error("Vector::Vector(gsl_matrix* vect) mismatch in sizes");
	}

}

Vector::Vector(Vector&& other) {

	vector = other.vector;
	orientation = other.orientation;
	length = -1;

	other.vector = nullptr;

}

Vector::Vector(const Vector& v) : orientation{v.orientation} {

	vector = create_vector();
	gsl_matrix_memcpy(vector,v.vector);
	length = -1;	

}

Vector Vector::operator*(const Matrix& rm) const {

	if(this->isvertical()){
		throw runtime_error("Vector::operator*(Matrix), wrong orientation of the vector");
	}

	Vector result;
	matmul(this->get_vector(),rm.get_matrix(),result.get_vector());
	return result;

}

Matrix Vector::operator*(const Vector& rv) const {

	if(this->ishorizontal() and rv.isvertical()){
		throw runtime_error("Vector::operator*(Vector), wrong orientation of the vectors");
	}

	Matrix result;
	matmul(this->get_vector(),rv.get_vector(),result.get_matrix());
	return result;

}


Vector& Vector::operator=(Vector&& other) {

	if(this == &other) return *this;

	gsl_matrix_free(vector);
	
	vector = other.vector;
	orientation = other.orientation;

	other.vector = nullptr;	

	length = -1;	

	return *this;

}

double Vector::v1() const {

	return gsl_matrix_get(vector,0,0);

} 

double Vector::v2() const {
	
	int i=1;
	int j=0;
	if(ishorizontal()) {
		i=0;
		j=1;
	}
	return gsl_matrix_get(vector,i,j);

}

double Vector::v3() const {

	int i=2;
	int j=0;
	if(ishorizontal()) {
		i=0;
		j=2;
	}
	return gsl_matrix_get(vector,i,j);

}

double Vector::get_length() {

	if(length == -1) {

		double a = v1();
		double b = v2();
		double c = v3();

		length = sqrt(a*a + b*b + c*c);

	}
	
	return length;

}

void Vector::normalize() {

	double d = get_length();

	if(d==0) {
		
		throw runtime_error("Vector::normalize(), length of the vecotr is equal to 0");

	}	

	double a = v1();
	double b = v2();
	double c = v3();

	if(ishorizontal()) {
		gsl_matrix_set(vector,0,0,a/d);
		gsl_matrix_set(vector,0,1,b/d);
		gsl_matrix_set(vector,0,2,c/d);
	}else {

		gsl_matrix_set(vector,0,0,a/d);
		gsl_matrix_set(vector,1,0,b/d);
		gsl_matrix_set(vector,2,0,c/d);

	}	

}

Vector Vector::give_transpose() const {

	Orientation ori = Orientation::horizontal;

	if(ishorizontal()) {

		ori = Orientation::vertical;

	}

	Vector result = Vector({v1(), v2(), v3()},ori);	

	return result;

}

Vector Vector::give_unit_vector(Orientation o) const {

	double a = v1();
	double b = v2();
	double c = v3();

	double d = sqrt(a*a + b*b + c*c);

	if(d==0) {

		throw runtime_error("Vector::give_unit_vector(), lenght of the vector is equal to 0");

	}

	Vector result = Vector({a/d,b/d,c/d},o);

	return result;

}

gsl_matrix* Vector::get_vector() const{

	return vector;

}

gsl_matrix* Vector::create_vector() {

	if(ishorizontal()) {
		length = -1;
		return gsl_matrix_alloc(1,3);
	}else if(isvertical()) {
		length = -1;
		return gsl_matrix_alloc(3,1);
	}else {
		throw runtime_error("Vector::create_vector(), not known orientation");
	}
}

bool Vector::ishorizontal() const {

	if(orientation==Orientation::horizontal) {
		return true;
	}else if(orientation==Orientation::vertical) {
		return false;
	}else {
		throw runtime_error("Vector::ishorizontal(), not known orientation");
	}

}

bool Vector::isvertical() const {

	if(orientation==Orientation::vertical) {
		return true;

	}else if(orientation==Orientation::horizontal) {
		return false;
	}else {
		throw runtime_error("Vector::isvertical(), not known orientation");
	}
}

void Vector::print() const {

	if(ishorizontal()) {
		for(int j=0;j<vector->size2;j++) {
			printf("%.2f ",gsl_matrix_get(vector,0,j));
        	}
		printf("\n");
	}else {
		for(int j=0;j<vector->size1;j++) {
			printf("%.2f\n",gsl_matrix_get(vector,j,0));
        	}
	}
}

Vector::~Vector() {

	if(vector!=NULL) {
		gsl_matrix_free(vector);
	}
}
