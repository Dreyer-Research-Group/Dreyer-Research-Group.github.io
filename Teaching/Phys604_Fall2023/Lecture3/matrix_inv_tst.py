import unittest

import numpy as np

def invert(A):
    ''' Function that finds the inverse of A '''
    
    if type(A) is not np.ndarray:
        raise Exception("Input must be an array")
    elif A.ndim != 2:
        raise Exception("Array must be 2D")
    elif A.shape[0] != A.shape[1]:
        raise Exception("Matrix must be square")

    try:

        Ainv=np.linalg.inv(A)
        
    except np.linalg.LinAlgError as err:
        if 'Singular matrix' in str(err):
            raise Exception("Matrix is singular")
        else:
            raise


    return Ainv
        
class TestInvMethods(unittest.TestCase):

    def test_inverse(self):

        # Some small number
        epsilon=1.0e-10
        
        # Make some test case:
        Amat=np.array([[1.0,2.0],[3.0,4.0]])
        avec=np.array([2.0,3.0])
        bvec=np.matmul(Amat,avec)

        diff_vec=np.linalg.norm(avec-np.matmul(invert(Amat),bvec))
        
        
        self.assertTrue(diff_vec<epsilon)

        
        
if __name__ == '__main__':
    unittest.main()
