
def refill_array(axis, dimension, new_array, old_array):

    if dimension == 1:
        N = old_array.size
        new_array[:N] = old_array
    elif dimension == 2:
        N = old_array.shape[axis]
        if axis == 0:
            for i in range(N):
                new_array[i,:] = old_array[i,:]
        elif axis == 1:
            for i in range(N):
                new_array[:,i] = old_array[:,i]
    elif dimension == 3:
        N = old_array.shape[axis]
        if axis == 0:
            for i in range(N):
                new_array[i,:,:] = old_array[i,:,:]
        elif axis == 1:
            for i in range(N):
                new_array[:,i,:] = old_array[:,i,:]
        elif axis == 2:
            for i in range(N):
                new_array[:,:,i] = old_array[:,:,i]
    elif dimension == 4:
        N = old_array.shape[axis]
        if axis == 0:
            for i in range(N):
                new_array[i,:,:,:] = old_array[i,:,:,:]
        elif axis == 1:
            for i in range(N):
                new_array[:,i,:,:] = old_array[:,i,:,:]
        elif axis == 2:
            for i in range(N):
                new_array[:,:,i,:] = old_array[:,:,i,:]
        elif axis == 3:
            for i in range(N):
                new_array[:,:,:,i] = old_array[:,:,:,i]