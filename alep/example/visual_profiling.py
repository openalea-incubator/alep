""" Graphical profiling of statement : 
    >>> statement = 'for i in range( 1000 ): svd( dot( standard_normal( ( 10, 200 ) ), 
                        standard_normal( ( 200, 10 ) ) ) )'
    >>>> profviz( statement )
    """

# def profviz(statement, filename = 'profile.pstats', sortkey = 'cumu' ):
    # from cProfile import run
    # from subprocess import call
    # run( statement, filename )
    # call( 'gprof2dot -f pstats output.pstats | dot -Tpng -o output.png', shell = True )
    
def profviz(statement, filename = 'profile', sortkey = 'cumu' ):
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from cProfile import run
    from subprocess import call
    from pstats import Stats
    run( statement, filename + '.pstats')
    Stats( filename + '.pstats').strip_dirs().sort_stats( sortkey ).print_stats(10)
    call( 'gprof2dot -f pstats ' + filename + '.pstats | dot -Tpng -o ' + filename + '.png', shell = True )
    img = mpimg.imread( filename + '.png'  )
    plt.imshow( img )
    
def profviz_notebook(statement, filename = 'profile.pstats', sortkey = 'cumu' ):
    from IPython.display import display, SVG
    from cProfile import run
    from subprocess import call
    from pstats import Stats
    run( statement, filename + '.pstats')
    Stats( filename + '.pstats').strip_dirs().sort_stats( sortkey ).print_stats(10)
    call( 'gprof2dot -f ' + filename + '.pstats | dot -Tsvg -o ' + filename + '.svg', shell = True )
    display( SVG( filename = filename + '.svg' ) )