from DXseq import *
import string

################################################################################
# imported variables
#dup1 = DXseq.dup1
#dup2 = DXseq.dup2
#hp1 = DXseq.hp1
#hp2 = DXseq.hp2
#tiletype = DXseq.tiletype
#tiletype2 = DXseq.tiletype2
#STACK_NUM = DXseq.STACK_NUM


# global variables
k=1
m=1
junction=[12,13,28,29]

# remove existing pdb files if they exist
if os.path.exists(DIR + r'/' + tiletype + '.pdb'):
    os.remove(tiletype + '.pdb')

pdbfile=open(tiletype + ".pdb","a")

################################################################################

def filewrite(coord,init_atomnum,fin_atomnum,strdnum,base,norm,i):

    global k,m
    
    if norm:
        for j in xrange(init_atomnum,fin_atomnum):

            pdbfile.write(('ATOM'+' '*2+'%5d'+' '*2+'%-3s'+' '*3+'%s'+' '+\
                               '%1s%4d'+' '*4+'%8.3f'*3+' '*2+\
                               '1.00'+' '*2+'0.00'+' '*11+'%s\n') %\
                              (k,atom_name[j],base,strdnum,m,coord[i][j][0],\
                                   coord[i][j][1],coord[i][j][2],atom_name[j][0]))

            k+=1
        m+=1

    else:
        if base=='A': base='T'
        elif base=='G': base='C'
        elif base=='C': base='G'
        elif base=='T': base=='A'

        for j in xrange(init_atomnum,fin_atomnum):
            
            pdbfile.write(('ATOM'+' '*2+'%5d'+' '*2+'%-3s'+' '*3+'%s'+' '+\
                               '%1s%4d'+' '*4+'%8.3f%8.3f%8.3f'+' '*2+\
                               '1.00'+' '*2+'0.00'+' '*11+'%s\n') %\
                              (k,atom_name[j],base,strdnum,m,coord[i][j][0],\
                                   coord[i][j][1],coord[i][j][2],atom_name[j][0]))
            
            k+=1
        m+=1

def DXoutput(coord1,coord2,tilenum,hpcoord1,hpcoord2,strdnum):

    # map letters to numbers
    for i in xrange(len(string.lowercase)):
        if strdnum==string.lowercase[i]:
            strdnum2=i+1

    if strdnum2%4==1:
        
        # strand 1
        # first sequence
        ibasenum=0
        fbasenum=junction[1]
        inc=1
        i=ibasenum
        
        for base in dup1[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord1,0,21,strdnum,base,True,i)
                i+=inc
            elif base=='G': 
                filewrite(coord1,41,63,strdnum,base,True,i)
                i+=inc
            elif base=='C': 
                filewrite(coord1,82,101,strdnum,base,True,i)
                i+=inc
            else: 
                filewrite(coord1,123,143,strdnum,base,True,i)
                i+=inc

        # strand 1        
        # second sequence
        ibasenum=junction[0]
        fbasenum=None
        inc=-1
        i=ibasenum

        if (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')\
                and tilenum == 'B':

            ibasenum=junction[0]
            fbasenum=4
            inc=-1
            i=ibasenum
            
        for base in dup2[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord2,21,41,strdnum,base,False,i)
                i+=inc
            elif base=='G': 
                filewrite(coord2,63,82,strdnum,base,False,i)
                i+=inc
            elif base=='C': 
                filewrite(coord2,101,123,strdnum,base,False,i)
                i+=inc
            else: 
                filewrite(coord2,143,164,strdnum,base,False,i)
                i+=inc

    elif strdnum2%4==2:

        # strand 2
        # first sequence        
        ibasenum=STACK_NUM-1
        fbasenum=junction[2]
        inc=-1
        i=ibasenum

        if tiletype == 'STLOX' or tiletype == 'STLXX':
            ibasenum=STACK_NUM-6
            fbasenum=junction[2]
            inc=-1
            i=ibasenum

        for base in dup1[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord1,21,41,strdnum,base,False,i)
                i+=inc
            elif base=='G': 
                filewrite(coord1,63,82,strdnum,base,False,i)
                i+=inc
            elif base=='C': 
                filewrite(coord1,101,123,strdnum,base,False,i)
                i+=inc
            else: 
                filewrite(coord1,143,164,strdnum,base,False,i)
                i+=inc

        # strand 2
        # second sequence
        ibasenum=junction[3]
        fbasenum=STACK_NUM
        inc=1
        i=ibasenum
        
        if (tiletype == 'STLOX' or tiletype == 'STLXX')\
                or ((tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')\
                        and tilenum == 'B'):
            ibasenum=junction[3]
            fbasenum=STACK_NUM-5
            inc=1
            i=ibasenum

        for base in dup2[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord2,0,21,strdnum,base,True,i)
                i+=inc
            elif base=='G': 
                filewrite(coord2,41,63,strdnum,base,True,i)
                i+=inc
            elif base=='C': 
                filewrite(coord2,82,101,strdnum,base,True,i)
                i+=inc
            else: 
                filewrite(coord2,123,143,strdnum,base,True,i)
                i+=inc
            
    elif strdnum2%4==3:

        # strand 3
        # first sequence
        ibasenum=5
        fbasenum=junction[3]
        inc=1
        i=ibasenum

        if (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')\
                and tilenum == 'B':
            ibasenum=0
            fbasenum=junction[3]-8
            inc=1
            i=ibasenum

        elif tiletype2 == 'DTLXO-2' and tilenum == 'B':
            ibasenum=5
            fbasenum=junction[3]-8
            inc=1
            i=ibasenum

        for base in dup2[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord2,0,21,strdnum,base,True,i)
                i+=inc
            elif base=='G': 
                filewrite(coord2,41,63,strdnum,base,True,i)
                i+=inc
            elif base=='C': 
                filewrite(coord2,82,101,strdnum,base,True,i)
                i+=inc
            else: 
                filewrite(coord2,123,143,strdnum,base,True,i)
                i+=inc


################################################################################
# hairpin1
        if (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXO-2'\
                or tiletype2 == 'DTLXX-2') and tilenum == 'B':
            # 5' -> 3' up the hairpin

            ibasenum=2
            fbasenum=12
            inc=1
            i=ibasenum
            
            for base in hp1[ibasenum:fbasenum:inc]:
                
                if base=='A':
                    filewrite(hpcoord1,0,21,strdnum,base,True,i)
                    i+=inc
                elif base=='G':
                    filewrite(hpcoord1,41,63,strdnum,base,True,i)
                    i+=inc
                elif base=='C':
                    filewrite(hpcoord1,82,101,strdnum,base,True,i)
                    i+=inc
                elif base=='T':
                    filewrite(hpcoord1,123,143,strdnum,base,True,i)
                    i+=inc

            # 5' -> 3' down the hairpin
            ibasenum=11
            fbasenum=None
            inc=-1
            i=ibasenum

            for base in hp1[ibasenum:fbasenum:inc]:
                
                if base=='A':
                    filewrite(hpcoord1,21,41,strdnum,base,False,i)
                    i+=inc
                elif base=='G':
                    filewrite(hpcoord1,63,82,strdnum,base,False,i)
                    i+=inc
                elif base=='C':
                    filewrite(hpcoord1,101,123,strdnum,base,False,i)
                    i+=inc
                elif base=='T':
                    if i==10 or i==11:
                        filewrite(hpcoord1,143,164,strdnum,base,True,i)
                        i+=inc
                    else: 
                        filewrite(hpcoord1,143,164,strdnum,base,False,i)
                        i+=inc

            # end of hairpin1
            # rest of strand
            ibasenum=junction[3]-8
            fbasenum=junction[3]
            inc=1
            i=ibasenum
            
            for base in dup2[ibasenum:fbasenum:inc]:
                
                if base=='A': 
                    filewrite(coord2,0,21,strdnum,base,True,i)
                    i+=inc
                elif base=='G': 
                    filewrite(coord2,41,63,strdnum,base,True,i)
                    i+=inc
                elif base=='C': 
                    filewrite(coord2,82,101,strdnum,base,True,i)
                    i+=inc
                else: 
                    filewrite(coord2,123,143,strdnum,base,True,i)
                    i+=inc

################################################################################

        # strand 3
        # second sequence
        ibasenum=junction[2]
        fbasenum=4
        inc=-1
        i=ibasenum
        
        for base in dup1[ibasenum:fbasenum:inc]:
            
            if base=='A': 
                filewrite(coord1,21,41,strdnum,base,False,i)
                i+=inc
            elif base=='G': 
                filewrite(coord1,63,82,strdnum,base,False,i)
                i+=inc
            elif base=='C': 
                filewrite(coord1,101,123,strdnum,base,False,i)
                i+=inc
            else: 
                filewrite(coord1,143,164,strdnum,base,False,i)
                i+=inc


    elif strdnum2%4==0:

        # strand 4
        # first sequence
        ibasenum=STACK_NUM-6
        fbasenum=junction[0]
        inc=-1
        i=ibasenum

        if tiletype == 'STLOX' or tiletype == 'STLXX':
            ibasenum=STACK_NUM-1
            fbasenum=junction[0]
            inc=-1
            i=ibasenum

        elif (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2')\
                and tilenum == 'B':
            ibasenum=STACK_NUM-1
            fbasenum=junction[0]+8
            inc=-1
            i=ibasenum

        elif tiletype2 == 'DTLXO-2' and tilenum == 'B':
            ibasenum=STACK_NUM-6
            fbasenum=junction[0]+8
            inc=-1
            i=ibasenum

        for base in dup2[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord2,21,41,strdnum,base,False,i)
                i+=inc
            elif base=='G': 
                filewrite(coord2,63,82,strdnum,base,False,i)
                i+=inc
            elif base=='C': 
                filewrite(coord2,101,123,strdnum,base,False,i)
                i+=inc
            else: 
                filewrite(coord2,143,164,strdnum,base,False,i)
                i+=inc

################################################################################
# hairpin2
        if (tiletype2 == 'DTLOX-2' or tiletype2 == 'DTLXX-2'\
                or tiletype2 == 'DTLXO-2') and tilenum == 'B':
            # 5' -> 3' up the hairpin
            ibasenum=2
            fbasenum=12
            inc=1
            i=ibasenum

            for base in hp2[ibasenum:fbasenum:inc]:
                
                if base=='A':
                    filewrite(hpcoord2,0,21,strdnum,base,True,i)
                    i+=inc
                elif base=='G':
                    filewrite(hpcoord2,41,63,strdnum,base,True,i)
                    i+=inc
                elif base=='C':
                    filewrite(hpcoord2,82,101,strdnum,base,True,i)
                    i+=inc
                elif base=='T':
                    filewrite(hpcoord2,123,143,strdnum,base,True,i)
                    i+=inc

            # 5' -> 3' down the hairpin
            ibasenum=11
            fbasenum=None
            inc=-1
            i=ibasenum

            for base in hp2[ibasenum:fbasenum:inc]:
                
                if base=='A':
                    filewrite(hpcoord2,21,41,strdnum,base,False,i)
                    i+=inc
                elif base=='G':
                    filewrite(hpcoord2,63,82,strdnum,base,False,i)
                    i+=inc
                elif base=='C':
                    filewrite(hpcoord2,101,123,strdnum,base,False,i)
                    i+=inc
                elif base=='T':
                    if i==10 or i==11:
                        filewrite(hpcoord2,143,164,strdnum,base,True,i)
                        i+=inc
                    else: 
                        filewrite(hpcoord2,143,164,strdnum,base,False,i)
                        i+=inc

            # end of hairpin2
            # rest of strand
            ibasenum=junction[0]+8
            fbasenum=junction[0]
            inc=-1
            i=ibasenum

            for base in dup2[ibasenum:fbasenum:inc]:            
                if base=='A':
                    filewrite(coord2,21,41,strdnum,base,False,i)
                    i+=inc
                elif base=='G':
                    filewrite(coord2,63,82,strdnum,base,False,i)
                    i+=inc
                elif base=='C':
                    filewrite(coord2,101,123,strdnum,base,False,i)
                    i+=inc
                else:
                    filewrite(coord2,143,164,strdnum,base,False,i)
                    i+=inc
################################################################################

        # strand 4
        # second sequence
        ibasenum=junction[1]
        fbasenum=STACK_NUM-5
        inc=1
        i=ibasenum
        
        if tiletype == 'STLOX' or tiletype == 'STLXX':
            ibasenum=junction[1]
            fbasenum=STACK_NUM
            inc=1
            i=ibasenum
        
        for base in dup1[ibasenum:fbasenum:inc]:

            if base=='A': 
                filewrite(coord1,0,21,strdnum,base,True,i)
                i+=inc
            elif base=='G': 
                filewrite(coord1,41,63,strdnum,base,True,i)
                i+=inc
            elif base=='C': 
                filewrite(coord1,82,101,strdnum,base,True,i)
                i+=inc
            else: 
                filewrite(coord1,123,143,strdnum,base,True,i)
                i+=inc


def Dvectorout(axis1,axis2,dup1,dup2):
    """returns the magnitude of the D-vector to 'Dvector.txt'"""

    # writeout Dvector
    Dvector=open("Dvector_" + tiletype  + ".txt",'w')
    Dcomp=open("Dcomp_" + tiletype  + ".txt",'w')

    origin = [0.,0.,0.]
    avec = [0. for x in xrange(3)]
    bvec = [0. for x in xrange(3)]
    abvec = [0. for x in xrange(3)]
    cvec = [0. for x in xrange(3)]
    dvec = [0. for x in xrange(3)]

    # define B-vector, difference of axes positions of the 
    # first base pairs between the two positions.
    for comp in xrange(3):
        bvec[comp] = axis2[0][comp] - axis1[0][comp]

    for base in xrange(1,STACK_NUM-5):

        for comp in xrange(3):

            avec[comp] = axis1[base][comp] - axis1[0][comp]
            abvec[comp] = avec[comp] + bvec[comp]
            cvec[comp] = axis2[base][comp] - axis1[0][comp]
            dvec[comp] = cvec[comp] - abvec[comp]

        dist = distance(dvec,origin)
        Dvector.write(('%3d' + ' '*5 + '%10f\n') % (base,dist))
        Dcomp.write(('%10f %10f %10f\n') % (dvec[0],dvec[1],dvec[2]))

    Dvector.close()
    Dcomp.close()
