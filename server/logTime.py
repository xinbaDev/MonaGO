def logTime(func):
    '''
    decorator

    Args:
        take the decorated func 
    
    Return:
        return the wrapped function

    '''
    def func_warpper(*arg, **kw):
        start_time = time.time()
        result = func(*arg, **kw)
        logger.info(func.__name__ + " lasts--- %s seconds ---" % (time.time() - start_time))
        return result
    return func_warpper
