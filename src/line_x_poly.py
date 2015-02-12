# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 21:37:13 2014
Modified from line_x_poly.pyx
@author: Weigang Wei
"""

import numpy as np
def intersection(x1, y1, x2, y2, x3, y3, x4, y4):
    ''' calcualte the intersection of two lines defined by piont1 ->x1, y1,, 
        point 2-> x2, y2, point 3-> x3, y3 and point4-> x4,y4. 
        p1 and p2 stand for line 1
        p3 and p4 stand for line 2
    '''
    try:
        x = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
        y = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
        return [x,y]
    except:
        # one point coincides each other 
        if (x1==x3 and y1==y3) or (x1==x4 and y1==y4):
            return [x1, y1]
        elif (x2==x3 and y2==y3) or (x2==x4 and y2==y4):
            return [x2,y2]
        else:
            return []
    

def is_interSegment(x, y, x1, y1, x2, y2, x3, y3, x4, y4):
    ''' Check wether point is the intersection point of line [lineStartPt, lineEndPt]
         and line [lineStartPt2, lineEndPt2]
        point = [x, y] 
        lineStartPt = [x, y]
        lineEndPt = [x, y]
        lineStartPt2 = [x, y]
        lineEndPt2 = [x, y]
    '''
    leftPtLi = min(x1, x2)
    rightPtLi = max(x1, x2)
    highPtLi = max(y1, y2)
    lowPtLi = min(y1, y2)
    leftPtL2 = min(x3, x4)
    rightPtL2 = max(x3, x4)
    highPtL2 =  max(y3, y4)
    lowPtL2 = min(y3, y4)

    xlimMin = max(leftPtLi, leftPtL2)
    xlimMax = min(rightPtLi, rightPtL2)
    ylimMin = max(lowPtLi, lowPtL2)
    ylimMax = min(highPtLi, highPtL2)
    if x>=xlimMin and x<=xlimMax and y>=ylimMin and y<=ylimMax:
        return True
    else:
        return False
        
def line_x_poly(lineStartPtx, lineStartPty, lineEndPtx, lineEndPty,\
              polygonx, polygony):
    ''' Find the intersection by a line and a polygon    
        WARNING: if the instersection point is coincide with the joint point 
        between two line segment, the resutls will return one point more than 
        once.         
        lineStartPt = [x, y]
        lineEndPt = [x, y]
        polygonx = [x1, x2, x3,..., xn, x1]
        polygony = [y1, y2, y3,...,yn, y1]        
    '''
    try:
        N = len(polygonx)-1
        intersecPt = []
        for n in xrange(N):
            x3, y3 = polygonx[n], polygony[n]
            x4, y4 = polygonx[n+1], polygony[n+1]
            ipt = intersection(lineStartPtx, lineStartPty, lineEndPtx, lineEndPty, x3, y3, x4, y4)        
            if ipt:
                isintsec = is_interSegment(ipt[0], ipt[1], lineStartPtx, lineStartPty, lineEndPtx, lineEndPty, x3, y3, x4, y4)
                if isintsec:
                    intersecPt.append(ipt)                     
        return intersecPt
    except:
        assert (len(polygonx)==len(polygony)), 'The size x-corrdinate must be the same with y-coordinate'
        assert (polygonx[0]==polygonx[-1] and polygony[0]==polygony[-1]), 'The start point of the vector must be coincide with the last point'
        raise
        return []

if __name__=='__main__':
    p1 = [0., 0.]
    p2 = [10., 10]
    polyx = np.array([1, 5, 5, 1, 1])
    polyy = np.array([1, 1, 10, 10, 1])
    intersecPt = line_x_poly(0., 0., 10.,10., polyx,polyy)
    print 'intersecPt'
    print intersecPt
        
        
