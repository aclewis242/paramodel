def procHex(*args):
    '''
    Converts the given values (0-255) into a hex color code.
    '''
    return '#'+''.join([f'0{hex(s)[2:]}'[-2:] for s in args])

def str2Color(s: str):
    '''
    Converts the given string into an arbitrary (but consistent) color code.
    '''
    tot_num = 3*2551*sum([ord(c)**5 for c in s])
    R = int(tot_num%255)
    G = int((tot_num/173)%255)
    B = int((tot_num*371)%255)
    if R + G + B > 2.5*255: [R, G, B] = [int(c/2) for c in [R, G, B]]
    return procHex(R, G, B)

def pop2Alpha(p: str):
    '''
    Returns an alpha value corresponding with the kind of population (S, I, R) submitted. (R is lighter.)
    '''
    p = p[0]
    if p == 'S' or p == 'I': return 1.
    else: return 0.7