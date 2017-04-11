def array_truncater(array):
    """In order to rebin by 5 pixels, count the length of the arrays 
    and delete data from extremeties such that the length is 5.
    
    Input: array - enter a 1-D array.
    Output: array - returns the truncated input array such that it is a multiple of 5."""
    
    count = len(array) % 5
    
    while(count != 0):
        array = np.delete(array, [0])
        count = count - 1
        
    return array