
import copy

def refill_array(axis, dimension, new_array, old_array):
    '''
    Used when iterating a dimension of an array to backfill the previously
    unelevated dimension. Mostly used within io.iterate_dimension()
    '''

    output = copy.deepcopy(new_array)

    if dimension == 1:
        N = old_array.size
        output[:N] = old_array
    elif dimension == 2:
        N = old_array.shape[axis]
        if axis == 0:
            for i in range(N):
                output[i,:] = old_array[i,:]
        elif axis == 1:
            for i in range(N):
                output[:,i] = old_array[:,i]
    elif dimension == 3:
        N = old_array.shape[axis]
        if axis == 0:
            for i in range(N):
                output[i,:,:] = old_array[i,:,:]
        elif axis == 1:
            for i in range(N):
                output[:,i,:] = old_array[:,i,:]
        elif axis == 2:
            for i in range(N):
                output[:,:,i] = old_array[:,:,i]
    elif dimension == 4:
        N = old_array.shape[axis]
        if axis == 0:
            for i in range(N):
                output[i,:,:,:] = old_array[i,:,:,:]
        elif axis == 1:
            for i in range(N):
                output[:,i,:,:] = old_array[:,i,:,:]
        elif axis == 2:
            for i in range(N):
                output[:,:,i,:] = old_array[:,:,i,:]
        elif axis == 3:
            for i in range(N):
                output[:,:,:,i] = old_array[:,:,:,i]

    return output