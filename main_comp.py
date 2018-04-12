import numpy as np 
import re, pprint, scipy, scipy.linalg
from numpy import linalg as LA
from numpy.linalg import inv
from scipy.linalg import hilbert, lu
from scipy.interpolate import lagrange
from fractions import Fraction


# get the shape of the data e.g. [10][5]
def getShape(data_str):
    count_semi = 0

    first_semi = data_str.find(';')                 # find first ';' in the line
    pos_semi = data_str[:first_semi].count(',')     # count amount of ',' from start of line to ';'
    count_semi += data_str.count(';')               # count amount of ';' in file
    
    return count_semi, pos_semi


# push contents of file to matrix
def pushToMatrix(data_str,i,j):
    current_digit = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", data_str)
    x = np.asarray(current_digit)
    ans = np.split(x,i)
    
    m = np.asmatrix(ans).astype(np.float64)
    return m


# convert matrix to comma-semi-colon-format
# print message to terminal & write output to file
def exportTomFormat(m):
    s = ""
    i,j = m.shape
    print "Matrix shape: [" + str(i) +" "+ str(j) +"]\n"
    for row in range(0,i):
        for col in range(0,j):
            s += str(m[row,col])
            if col == j-1:
                s
            else:
                s += ","

        s += ";"

    # write output to a file (will create file if not existing)
    # NOTE: does NOT append to file - will overwrite contents
    f = open('output_file.txt', 'w')
    f.write(s)  # python will convert \n to os.linesep
    f.close()

    print "- Written to 'output.txt' -\n\n**************************************\n"


# # sum contents of matrix
# def addMatrixElement(m):
#     mat = m.flatten().sum()
#     return mat

# sum each row
def sumOfRow(m):
    m = np.asmatrix(m)
    i,j = m.shape
    arr = []

    for row in m:
        r = 0
        r = np.sum(row)
        arr.append(r)

    return arr

# get euclidean norm of matrix
def vectorNorm(m):
    mat = m.flatten()
    norm_matrix = LA.norm(mat)                  # Euclidean norm
    return norm_matrix


# using numpy library for cond and infinity norm
def condNumberNumpy(m):
    inf_norm_matrix = LA.cond(m, np.inf)        # infinity norm inf = max(sum(abs(x), axis=1))
    return inf_norm_matrix


# own cond number fuction
# needs to be passed a matrix
def condNumber(m):
    matrix_total = 0
    inverse_matrix_total = 0


    # matrix
    for i in m:
        print "cond of row: " + str(abs(i).sum())
        if matrix_total < abs(i).sum():
            matrix_total = abs(i).sum()                      
    print "The largest: [" + str(matrix_total) + "]\n"

    # inverse of matrix
    inv_m = inv(m)
    for i in inv_m:
        print "cond of row: " + str(abs(i).sum())
        if inverse_matrix_total < abs(i).sum():
            inverse_matrix_total = abs(i).sum()
    print "The largest: [" + str(inverse_matrix_total) + "]\n"
    
    return matrix_total * inverse_matrix_total


# convert matrix elements to fractions
def convertToFractionMatrix(decimal_matrix):
    two_d_fract_array = []
    x,y = decimal_matrix.shape
    for h in decimal_matrix:
        for i in h:
            two_d_fract_array.append(Fraction(str(i)))
    fract_matrix = np.asmatrix(np.split(np.asarray(two_d_fract_array),x)) # reshape the matrix
    return fract_matrix


# create Hilbert Matrix return matrix -> fractions
# Not great - better to use convertToFractionMatrix() after calculation
def createHilbertFractionMatrix(sq_size):
    fract_hill =[]
    hill = hilbert(sq_size)                    # for decimals just use this line (from library)
    for h in hill:
        for i in h:
            fract_hill.append(Fraction(str(i)))

    reshaped_array = np.split(np.asarray(fract_hill),sq_size)
    print "The Hilbert Matrix shape: " + str(np.asarray(reshaped_array).shape) + "\n"
    return np.asmatrix(reshaped_array)

# def gausselim(A,B):
#     """
#     Solve Ax = B using Gaussian elimination and LU decomposition.
#     A = LU   decompose A into lower and upper triangular matrices
#     LUx = B  substitute into original equation for A
#     Let y = Ux and solve:
#     Ly = B --> y = (L^-1)B  solve for y using "forward" substitution
#     Ux = y --> x = (U^-1)y  solve for x using "backward" substitution
#     :param A: coefficients in Ax = B
#     :type A: numpy.ndarray of size (m, n)
#     :param B: dependent variable in Ax = B
#     :type B: numpy.ndarray of size (m, 1)
#     """
#     # LU decomposition with pivot
#     pl, u = lu(A, permute_l=True)
#     # forward substitution to solve for Ly = B
#     y = np.zeros(B.size)
#     for m, b in enumerate(B.flatten()):
#         y[m] = b
#         # skip for loop if m == 0
#         if m:
#             for n in xrange(m):
#                 y[m] -= y[n] * pl[m,n]
#         y[m] /= pl[m, m]

#     # backward substitution to solve for y = Ux
#     x = np.zeros(B.size)
#     lastidx = B.size - 1  # last index
#     for midx in xrange(B.size):
#         m = B.size - 1 - midx  # backwards index
#         x[m] = y[m]
#         if midx:
#             for nidx in xrange(midx):
#                 n = B.size - 1  - nidx
#                 x[m] -= x[n] * u[m,n]
#         x[m] /= u[m, m]
#     return x

# get permutation from LU
def getPermutation(P):
    perm = []
    for i in P:
        for j in i:
            if j == 1:
                 perm.append(i.tolist().index(j)+1)

    return perm

# three column vectors from tridiagonal matrix
def getTridiagonal(m):
    print m
    x,y = m.shape
    vector_one = []
    vector_two = []
    vector_three = []
    i = 0

    while i < y and i < x:
        
        vector_two.append(m[i][i])
        if i < x-1:
            vector_three.append(m[i+1][i])
        if i < y-1:
            vector_one.append(m[i][i+1])
        i+=1

    return np.asmatrix(vector_one).T, np.asmatrix(vector_two).T, np.asmatrix(vector_three).T

def countNumLen(n):
    count=0
    while(n>0):
        count=count+1
        n=n//10



if __name__ == '__main__':
    f = open('read.txt', 'r').read().strip()    # save file to string

    # print "**************************************\n"

    # x,y = getShape(f)
    # print "Data shape: ["+ str(x) +"]["+ str(y) +"]\n"
    # print "**************************************\n"

    # read_matrix = pushToMatrix(f,x,y)
    # print "The data provided as a matrix: \n"+str(read_matrix)+"\n"
    # print "**************************************\n"

    # sum_matrix = addMatrixElement(read_matrix)
    # print "Elements added together: ["+str(sum_matrix)+"]\n"
    # print "**************************************\n"

    # norm_matrix = vectorNorm(read_matrix)
    # print "Vector Norm: ["+str(norm_matrix)+"]\n"
    # print "**************************************\n"

    # cond_norm_matrix = condNumberNumpy(read_matrix)
    # print "Condition number: ["+str(cond_norm_matrix)+"]\n"
    # print "**************************************\n"

    # print "Condition number extra input test: [[.01, 98796],[.0001, 7896]]\n"
    # print "Cond = " + str(condNumber(np.asmatrix([[.01, 98796],[.0001, 7896]]))) + "\n"
    # print "**************************************\n"

    # own_cond_norm_matrix = condNumber(read_matrix)
    # print "Own condition number: ["+str(own_cond_norm_matrix)+"]\n"
    # print "**************************************\n"

    # hillbert_matrix = createHilbertFractionMatrix(3)
    # print "Hilbert Matrix fractions: \n" + str(hillbert_matrix) + "\n"
    # print "**************************************\n"

    # sum_hilbert = addMatrixElement(hillbert_matrix)
    # print "Sum of Hilbert Matrix fractions: [" + str(sum_hilbert) + "]\n"
    # print "**************************************\n"

    # print "Multiply Hilbert by inverse\n"
    # print np.round(np.dot(hilbert(4),inv(hilbert(4))))
    # print
    # print "**************************************\n"

    # print "LU Decomposition: \n"
    # LU_array = np.asarray([[1,2,3],[4,5,6],[7,8,9]])
    # b = [[6,15,24]]

    # # P = permutation matrix
    # P, L, U = scipy.linalg.lu(LU_array, permute_l=False, overwrite_a=False, check_finite=True)

    # print "Permutation Matrix: \n" + str(P) +"\n"
    # print "Permutation : " + str(getPermutation(P)) + "\n"

    # # L = lower triangle
    # print "Lower triangle: \n" + str(L) +"\n"

    # # U = upper triangle
    # print "Upper triangle: \n" + str(U) +"\n"

    # bb = P*b
    # print "bb: \n" + str(np.asmatrix(bb))

    # y = np.linalg.lstsq(L,bb)
    # print "y: \n" + str(np.asmatrix(y))
    # print "**************************************\n"

    # top_diagonal, main_diagonal, bottom_diagonal = getTridiagonal(LU_array)
    # print
    # print "Tridiagonal of Matrix - three column vectors: \n\n" + str(top_diagonal) + "\n\n" + str(main_diagonal) + "\n\n" + str(bottom_diagonal) +"\n"
    # print "**************************************\n"

    # print main_diagonal.shape
    # print bottom_diagonal.shape
    # # print np.dot(B,c)

    # print
    # print main_diagonal.shape
    # print bottom_diagonal.shape
    # print np.cross(main_diagonal.T, bottom_diagonal.T)
    # print 
    # print "**************************************\n"

    # exportTomFormat(LU_array)

    # print "**************************************\n"
    # print "Cholesky: \n"

    # chol_matrix = [[4, 2, 1, 1],[2, 5, 0, 3],[1, 0, 7, 2],[1, 3, 2, 8]]
    # print np.matrix(chol_matrix)
    # print
    # # b = [[8,10,10,14]]
    # # print b
    # # print
    # c = np.linalg.cholesky(chol_matrix)
    # print c

    # print "\n**************************************\n"
    # print "Testing Hilbert: \n"
    # print condNumberNumpy(hilbert(10))
    # print LA.norm(hilbert(10), np.inf)
    # print "\n**************************************\n"
    # print "LaGrange Polynomial Interpolation: \n"
    # x = np.array(read_matrix.flatten())
    # print x
    # y = np.array(read_matrix.flatten())
    # poly = lagrange(x, y)
    # print poly
    # print "\n**************************************\n"
    # LU = np.asarray([[1,2,3],[4,5,6],[7,8,9]])
    # print condNumber(hilbert(7))

    for i in range(1,10):
        print "Hilbert = " + str(i) 
        A = hilbert(i)
        print "the hilbert matrix: \n"+str(A)+"\n"
        b = sumOfRow(hilbert(i))
        print "b: "
        print b
        x = np.linalg.lstsq(A,b)
        print "AX = b: "
        print x
        print




