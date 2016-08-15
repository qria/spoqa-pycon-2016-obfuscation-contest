import sys
DEBUG = False

def log(*args, **kwargs):
    """ Simple debugging. Use this instead of `print`, please. """
    if not DEBUG:
        return
    print(*args, **kwargs)

class TimeoutError(Exception):
    pass

def calc_power_series(coefficients, x):
    """ Calculate value of power series with given coefficients.
        Supports approximating infinite power series with
        Domb-Sykes method.
        (For details, see https://en.wikipedia.org/wiki/Radius_of_convergence#Domb.E2.80.93Sykes_plot) 
        
        If `coefficients` is a list or tuple or any iterable that implements `__len__` or `__length_hint__`, 
            this assumes finiteness and returns the exact value of power series.
        If not, assumes infiniteness and returns estimation.


        result = c_0 + c_1 * x + c_2 + x ** 2 + ...
        coefficients: (aliased as `coeffs` in code) 
            list (or generator) of numbers

        """

    coeffs = coefficients # alias because too long
    coeffs_cache = [] # cache to avoid calling generator twice
    # There's probably better ways of caching 

    # Try to guess if coefficient is finite 
    # Supports list, tuple, etc.
    is_finite = hasattr(coeffs, '__len__') or hasattr(coeffs, '__length_hint__')

    # Validation for inifinite power series
    # Estimate radius of convergence and test if x is in range of (-1/r, 1/r)  
    r = float('inf') # init r
    err = "x is (-r, r) or coefficients is list or finite"
    if not is_finite:
        # Domb-Sykes radius of convergence estimation
        radius_estimation = 10
        for _ in range(radius_estimation):
            coeffs_cache.append(next(coeffs))

        X = [] # inverse of n
        Y = [] # c_n / c_(n-1)
        for n in range(1, len(coeffs_cache)):
            X.append(1/(n+1))
            Y.append(coeffs_cache[n] / coeffs_cache[n-1])

        line = fit(X, Y)
        r = 1 / line(0) # estimated radius of convergence
        if r > 0 and not -r < x < r:
            raise ValueError(err)

    if is_finite: 
        # Finite power series is trivial to calculate
        finite = coeffs[0]
        for i, c in enumerate(coeffs[1:]):
            finite += c * x ** (i+1)
        result = finite
    else:
        infinite = 0
        threshold = 1000 # If it takes too long, stop
        err = 0.0000001 # Error boundary
        i = 0 
        for c in coeffs_cache:
            last_infinite = infinite
            infinite += c * x ** i
            i += 1

        try:
            while abs(infinite - last_infinite) > err:
                last_infinite = infinite
                c = next(coeffs)
                infinite += c * x ** i
                i += 1
                if(i > threshold):
                    raise TimeoutError() 
        except TimeoutError:
            log('TimeoutError Raised!')
            del infinite # in case I accidently use this
        else:
            result = infinite


    # For debugging
    if not isinstance(coeffs, (list, tuple)) and "infinite" not in locals():
        # Oh no, too much iteration in infinite 
        suggested_fix_template = "10 * {error}"
        log("Too much iteration! You might want to adjust error percision with formula %s" % suggested_fix_template)
        # Render template
        suggested_error_precision = eval(suggested_fix_template.format(error=err))
        log("which is %s" % suggested_error_precision)

        #sanity_check_template = "100 * {suggested_fix}"
        sanity_check_template = "1 * {suggested_fix}" # 100 was too big. Needs more testing
        log("If %s is bigger than 1, there's something wrong." % sanity_check_template)
        log("Note that Domb-Sykes method assumes the sequence to have same or alternating signs")
        sanity_check = eval(sanity_check_template.format(suggested_fix=suggested_error_precision)) # Render template
        print(sanity_check)

    return result


def fit(X, Y):
    """ Simple linear regression """
    mean_x = sum(X) / len(X)
    mean_y = sum(Y) / len(Y)

    b = sum((x - mean_x) * (y - mean_y) for (x, y) in zip(X, Y)) / sum((x-mean_x)**2 for x in X)
    a = mean_y - b * mean_x
    return lambda x: a + b*x


if __name__ == '__main__':
    """ Run this to test this script """

    # Test fit
    assert(fit([1, 2, 3, 4], [5, 6, 7, 8])(9) == 13) # y = x + 4
    assert(fit([1, 3, 2, 4], [5, 6, 7, 8])(9) == 11.7) # y = 0.8x + 4.5

    # Test calc_power_series
    ## finite series
    assert(calc_power_series([1, 2, 3, 4], 5) == 586) # list
    assert(calc_power_series((1, 2, 3, 4), -1) == -2) # tuple also works
    assert(calc_power_series([1, 2, 3, 4], 0.1) == 1.234) # float works 

    ## infinite series
    def yield_infinite_ones():
        # cf: 1+ x+ x**2 + x**3 + ... == 1 / (1-x) 
        while True:
            yield 1

    epsilon = 0.001 # precision
    assert(abs(calc_power_series(yield_infinite_ones(), 0.9) - 10) < epsilon) # 1 / (1-0.9) == 10

    ## x in invalid range 
    try:
        calc_power_series(yield_infinite_ones(), 1.1) # 1 / (1 - 1.1) < 0 !!!
    except ValueError:
        pass
    else:
        print('x in invalid range didn\'t raise an exception!')
        assert(False)


    if len(sys.argv) > 0:
        calc_power_series(sys.argv[1], 1)