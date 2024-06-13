#include <iostream>
#include <utility>
#include <memory>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "../includes/vpsc_classes.h"
#include "../includes/vpsc_util.h"

using namespace std;

Matrix::Matrix() {

	matrix = gsl_matrix_alloc(3,3);
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
			gsl_matrix_set(matrix,i,j,0);
		}
	}

}

Matrix::Matrix(initializer_list<double> arr) {

	matrix = gsl_matrix_alloc(3,3);
	
	initializer_list<double>::iterator itr = arr.begin();
	
	try{
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++) {
				gsl_matrix_set(matrix,i,j,*itr);
				itr++;
			}
		}
	}catch(runtime_error e) {
		cout<<e.what()<<endl;
	}

	if(itr!=arr.end()) {

		throw runtime_error("Matrix::Matrix(initializer_list): to many elements in array");

	}
}

Matrix::Matrix(gsl_matrix *mat) {

	matrix = gsl_matrix_alloc(3,3);
	gsl_matrix_memcpy(matrix,mat);

}

Matrix::Matrix(const Matrix& mat) {
	
	matrix = gsl_matrix_alloc(3,3);
	gsl_matrix_memcpy(matrix,mat.get_matrix());

}

Matrix::Matrix(Matrix&& other) {
	
	matrix = other.matrix;
	other.matrix = nullptr;

}

Matrix Matrix::operator*(const Matrix& rm) const {

	Matrix result;
	matmul(this->get_matrix(),rm.get_matrix(),result.get_matrix());
	return result;

}

Vector Matrix::operator*(const Vector& rv) const {

	if(rv.ishorizontal()) {
		throw runtime_error("Matrix::operator*(Vector), vector must be vertical");

	}

	Vector result(Vector::Orientation::vertical);
	matmul(this->get_matrix(),rv.get_vector(),result.get_vector());
	return result;

}

Matrix operator*(double d, const Matrix& mr) {


	Matrix result;

	for(int i=0;i<mr.get_matrix()->size1;i++) {
		for(int j=0;j<mr.get_matrix()->size2;j++) {
			gsl_matrix_set(result.get_matrix(),i,j,d*mr.get(i,j));		
		}
	}

	return result;

}

Matrix Matrix::operator+(const Matrix& mr) const {

	Matrix result;
	for(int i=0;i<matrix->size1;i++) {
		for(int j=0;j<matrix->size2;j++) {
			gsl_matrix_set(result.get_matrix(),i,j,this->get(i,j)+mr.get(i,j));
		}
	}

	return result;

}

Matrix Matrix::operator-(const Matrix& mr) const {

	Matrix result;
	for(int i=0;i<matrix->size1;i++) {
		for(int j=0;j<matrix->size2;j++) {
			gsl_matrix_set(result.get_matrix(),i,j,this->get(i,j)-mr.get(i,j));
		}
	}
	
	return result;

}

Matrix& Matrix::operator=(Matrix&& other) {

	if(this == &other) return *this;

	matrix = other.matrix;
	other.matrix = nullptr;	

	return *this;

}

double Matrix::get(size_t i, size_t j) const {

	if(i<0 or i>2 or j<0 or j>2) {
	
		throw runtime_error("Matrix::get(), wrong parameter for get element from matrix");

	}

	return gsl_matrix_get(matrix,i,j);

}

double Matrix::m11() const {

	return gsl_matrix_get(matrix,0,0);

} 

double Matrix::m12() const {

	return gsl_matrix_get(matrix,0,1);

}

double Matrix::m13() const {

	return gsl_matrix_get(matrix,0,2);

}

double Matrix::m21() const {

	return gsl_matrix_get(matrix,1,0);

}

double Matrix::m22() const {

	return gsl_matrix_get(matrix,1,1);

}

double Matrix::m23() const {

	return gsl_matrix_get(matrix,1,2);

}

double Matrix::m31() const {

	return gsl_matrix_get(matrix,2,0);

}

double Matrix::m32() const {

	return gsl_matrix_get(matrix,2,1);

}

double Matrix::m33() const {

	return gsl_matrix_get(matrix,2,2);

}

Matrix Matrix::give_inverse() const {

	int s; 

	Matrix invmat;
        gsl_permutation* p = gsl_permutation_alloc(3); 
 
        gsl_linalg_LU_decomp (this->matrix, p, &s);     
        gsl_linalg_LU_invert (this->matrix, p, invmat.get_matrix());

        gsl_permutation_free(p);

        return invmat;

}

Matrix Matrix::give_transpose() const {

	Matrix result;
	for(int i=0;i<matrix->size1;i++) {
		for(int j=0;j<matrix->size2;j++) {
			gsl_matrix_set(result.get_matrix(),i,j,get(j,i));

		}
	}

	return result;

}

gsl_matrix* Matrix::get_matrix() {

	return matrix;

}

gsl_matrix* Matrix::get_matrix() const {

	return matrix;

}


void Matrix::print() const {

        for(int i=0;i<matrix->size1;i++) {
                for(int j=0;j<matrix->size2;j++) {
                        printf("%.4f ",gsl_matrix_get(matrix,i,j));
                }
                printf("\n");
        }

}

Matrix::~Matrix() {

	if(matrix!=NULL) {
		gsl_matrix_free(matrix);
	}
}
