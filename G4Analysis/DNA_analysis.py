#1/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created at 2012.03.06
'''

import Simple_atom
import time as Time
import sys
import os
import usage
import MDAnalysis
import DNA_param
from math import cos, sin, sqrt
import math
import numpy
from numpy import matrix
from numpy import dot

import DNA_matrix
from Coor import atomlib

def Get_parm_fromTRJ(traj_file, coor_file, base_list_1, base_list_2, output_name,CALCU="step",skip=1, dt=1,begin=0,end=-1):
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    if len(base_list_1)==len(base_list_2):
        LIST_NUM=len(base_list_1)
    else:
        print "ERROR: The length of the base list not match."
        return -1
    if len(output_name)!=LIST_NUM:
        print "ERROR: The number of the output file not match the size of the base list."

    Atom_list=Simple_atom.Get_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for i in range(LIST_NUM):

        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp = open(output_name[i], 'w')
        fp.write("#Group 1: ")
        for j in base_list_1[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#Group 2: ")
        for j in base_list_2[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        if CALCU=="step":
            fp.write("#  shift       slide        rise        tilt        roll       twist\n")
        else:
            fp.write("#  shear     stretch     stagger      buckle   propeller     opening\n")
        fp.close()

#        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])

        temp_list=list()
        for m in base_list_1[i]:
            temp_list.append(residue_list[m-1])
        base_name_list_1.append(temp_list)
        # base_name_list_1.append(residue_list[m-1])

        temp_list=list()
        for m in base_list_2[i]:
            temp_list.append(residue_list[m-1])
        base_name_list_2.append(temp_list)


        base_atom_list_1.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1[i]])
        base_atom_list_2.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2[i]])
        # print base_atom_list_1


    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if dt > 0.0:
            if time >= float(begin):
                continue
            if time > float(end) and end !=-1:
                break

        if ts.frame % skip == 0 :
            for i in range(LIST_NUM):
                r1=[]
                '''the group 1 rotate list'''
                r2=[]
                '''the group 2 rotate list'''
                c1=[]
                '''the group 1 coordinate list'''
                c2=[]
                '''the group 2 coordinate list'''
                for m in range(len(base_name_list_1[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_1[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[i][m][0])
                    #base_name_list_1[index of the groups][index of the base of group 1][base_name,base_serial]
                    r1.append(result[0])
                    c1.append(result[1])

                for m in range(len(base_name_list_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
                    r2.append(result[0])
                    c2.append(result[1])

                fp = open(output_name[i], 'a')

                if CALCU=="pair":
                    a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))
                else:
                    middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
                    middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
                    a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
                    fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))

                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

                fp.close()

    print "The DNA helical analysis finished"
    print "The result are in file: %s" %output_name

def Get_para_fromTOP( coor_file, base_list_1, base_list_2,CALCU="step",PRINT=True):
    '''
    get pair parameters from coor_file
    '''

    atom_list=Simple_atom.Get_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(atom_list)
    # print residue_list
    # atom_list=Simple_atom.Get_atom_list(coor_file)

    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for m in base_list_1:
        # for n in residue_list:
        #     if n[1]==m:
        base_name_list_1.append(residue_list[m-1])

    # print base_name_list_1

    for m in base_list_2:
        base_name_list_2.append(residue_list[m-1])

    base_atom_list_1=[DNA_matrix.Get_baseID_list(atom_list,j) for j in base_list_1]
    base_atom_list_2=[DNA_matrix.Get_baseID_list(atom_list,j) for j in base_list_2]

    r1=[]
    '''the group 1 rotate list'''
    r2=[]
    '''the group 2 rotate list'''
    c1=[]
    '''the group 1 coordinate list'''
    c2=[]
    '''the group 2 coordinate list'''
    for m in range(len(base_name_list_1)):
        temp_list = [ [atom_list[x-1].atom_coor_x*10, atom_list[x-1].atom_coor_y*10,atom_list[x-1].atom_coor_z*10] \
                for x in base_atom_list_1[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[m][0])
        r1.append(result[0])
        c1.append(result[1])

    for m in range(len(base_name_list_2)):
        temp_list = [ [atom_list[x-1].atom_coor_x*10, atom_list[x-1].atom_coor_y*10,atom_list[x-1].atom_coor_z*10] \
                for x in base_atom_list_2[m] ]
#        print temp_list
        result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[m][0])
        r2.append(result[0])
        c2.append(result[1])

    if CALCU=="pair":
        a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
        if PRINT:
            print "   shear     stretch     stagger      buckle   propeller     opening"
            print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])
        return a
    else:
        middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
        middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])
        a=DNA_param.base_step_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
        if PRINT:
            print "   shift       slide        rise        tilt        roll       twist"
            print "%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f" %(a[0],a[1],a[2],a[3],a[4],a[5])
        return a


def Get_Dihedral_fromTOP( coor_file, base_list_1,PRINT=False ):
    '''
    calculate the dihedral from the top file.

    '''
    Atom_list=Simple_atom.Get_atom_list(coor_file)
    reside_list=Simple_atom.Get_Residue_list(Atom_list)

    if PRINT == True:
        print "%11s%10s%10s%10s%10s%10s%10s%10s%8s%8s%8s%8s%8s%8s%8s"\
                %("baseID","alpha","beta","gamma","delta","epslon","zeta","chi",\
                "alpha","beta","gamma","delta","epslon","zeta","chi")

    for m in base_list_1:
        resu=DNA_param.Get_Dihedral(Atom_list,m)   

        if PRINT == True:
            print "%8d%3s" %(reside_list[m-1][1],reside_list[m-1][0]),
            for i in range(7):
                if resu[i]=="-":
                    print "%9s" %("-"*4), 
                else:
                    print "%9.2f" %(resu[i]),
            for i in range(7):
                if resu[i]=="-":
                    print "%7s" %("-"*4), 
                else:
                    print "%7s" %(resu[i+7]),

            print "" 

#a        return [alpha,beta,gamma,delta,epslon,zeta,chi]


def Get_Dihedral_fromTRJ(traj_file, coor_file, base_list, output_name,skip=1, dt=1,begin=0,end=-1):
    '''
    note: if the atom index in coor file not start with 1, there will be an error.
    '''
    START_TIME=Time.time()

    Atom_list=Simple_atom.Get_atom_list(coor_file)


    for file_name in output_name:

        if os.path.isfile(file_name):
            print "backup %s to %s" %(file_name,"#"+file_name+"#")
            try:
                os.rename(file_name,"#"+file_name+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %file_name

        fp = open(file_name, 'w')
        fp.write("#Residue Index: ")
        # print base_list[output_name.index(file_name)]
        fp.write("%d\n" %base_list[output_name.index(file_name)])
        fp.write("#skip:%d\n" %skip)

        fp.write("%12s%10s%10s%10s%10s%10s%10s%10s%8s%8s%8s%8s%8s%8s%8s\n"\
                %("Time(ns)","alpha","beta","gamma","delta","epslon","zeta","chi",\
                "alpha","beta","gamma","delta","epslon","zeta","chi"))
        fp.close()

    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:
        time=float((ts.frame-1)*dt)
        if dt > 0.0:
            if time < float(begin):
                continue
            if time > float(end) and end !=-1:
                break

        if ts.frame % skip == 0 :
            for atom in Atom_list:
                atom.atom_coor_x=ts._x[atom.atom_id-1]
                atom.atom_coor_y=ts._y[atom.atom_id-1]
                atom.atom_coor_z=ts._z[atom.atom_id-1]


            for i in range(len(base_list)):
                fp = open(output_name[i], 'a')
                resu=DNA_param.Get_Dihedral(Atom_list,base_list[i])
                fp.write("%12.4f" %(time/1000))
                for j in range(7):
                    if resu[j]=="-":
                        fp.write("%10s" %("-"*4)) 
                    else:
                        fp.write("%10.2f" %(resu[j]))

                for j in range(7):
                    if resu[j]=="-":
                        fp.write("%8s" %("-"*4)) 
                    else:
                        fp.write("%8s" %(resu[j+7]))

                fp.write("\n") 


                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    if time < 1000:
                        usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))
                    elif time > 1000 and ts.frame %100 ==0 :
                        usage.echo("  analysis frame %6d, time %8.2f ns, time used %8.2f s\r" %(ts.frame, time/1000,NOW_TIME-START_TIME))

                fp.close()

    print "\nThe DNA helical analysis finished"
    print "The result are in file: %s" %output_name


def Get_pls_fromTRJ(traj_file, coor_file, base_list_1, base_list_2, output_name,skip=1, dt=1,begin=0,end=-1):
#    if coor_file.endswith('.top'):
#        Atom_list=amber_top
    START_TIME=Time.time()
    if len(base_list_1)==len(base_list_2):
        LIST_NUM=len(base_list_1)
    else:
        print "ERROR: The length of the base list not match."
        return -1
    if len(output_name)!=LIST_NUM:
        print "ERROR: The number of the output file not match the size of the base list."

    Atom_list=Simple_atom.Get_atom_list(coor_file)
    residue_list=Simple_atom.Get_Residue_list(Atom_list)
    chain=[]
    for ii,reside in enumerate(residue_list):
        if reside[0] in atomlib.RESIDUE_NAME_LIST:
            chain.append([ii+1,reside])  
        else:
            pass
    for ca in chain:
        print "%4d\t (  %6d  %8s )" % (ca[0], ca[1][1],ca[1][0])

    if (len(chain) % 2) == 0:
        pass
    else:
        print "-the base number of nucleic acid segment is not even, so..."
        print "-please check it!"


    #####
    # construct the lists for bp of nucleic acids 
    chain_bpList_1=list()
    chain_bpList_2=list()
    chain_bpAtomList_1=list()
    chain_bpAtomList_2=list()
    ##
    # chain_bpList_1 -> save the first bp of a neigbhor base pair step: 
    #    format: [ [[residue_name, residue_serial], [residue_name, residue_serial]], [[], []], ...]
    # chain_bpList_2 -> save the second bp of a neigbhor base pair step; 
    #    format: [ [[[residue_name, residue_serial], [residue_name, residue_serial]], [[], []], ... ]
    ##
    num_nuc=len(chain)
    for ii in range(num_nuc/2-1):
        # if (ii >= 0) and (ii < num_nuc/2-1):
        chain_bpList_1.append([ chain[ii][1], chain[num_nuc-ii-1][1] ])
        chain_bpList_2.append([ chain[ii+1][1], chain[num_nuc-ii-2][1] ])
        
    if len(chain_bpList_1) == len(chain_bpList_2) == num_nuc/2-1:
        chainLIST_NUM=len(chain_bpList_1)
        pass
    else :
        print "- number of base pairs in chain_bpList_1: %d" %(len(chain_bpList_1))
        print "- number of base pairs in chain_bpList_2: %d" %(len(chain_bpList_2))
        print "- total base pairs in system: %d" %(num_nuc/2)
        sys.exit()

    for ii, bp_1 in enumerate(chain_bpList_1):
        print ii, bp_1, chain_bpList_2[ii]

    #build the atom list for each base pair
    for i in range(chainLIST_NUM):
        chain_bpAtomList_1.append([DNA_matrix.Get_baseID_list(Atom_list,j[1]) for j in chain_bpList_1[i]])
        chain_bpAtomList_2.append([DNA_matrix.Get_baseID_list(Atom_list,j[1]) for j in chain_bpList_2[i]])
    # print chain_bpAtomList_1[0]
    # print chain_bpAtomList_2[0]
    # sys.exit()

    #####
    #build the list (base_list/base_atom_list) for the base list set in the parameter input file
    #base_name_list -> the residue list for base pairs which is gonna to be calculated
    # format: base_name_list_1: [ [[residue1_name, residue1_serial], [residue2_name, residue2_serial]], [[], []], ... ]
    # format: base_name_list_2: [ [[residue3_name, residue3_serial], [residue4_name, residue4_serial]], [[], []], ... ]
    # residue 1-2 -> bp-1
    # residue 3-4 -> bp-2
    ##
    #base_atom_list -> the 9/6 atoms (used in least square fitting) for each base
    # format: base_atom_list_1: [ [[1.1, 1.2, ..., 2.6/9 ], [2.1, 2.2, ..., 2.6/9 ]], [], ...]
    # format: base_atom_list_2: [ [[3.1, 3.2, ..., 2.6/9 ], [4.1, 4.2, ..., 2.6/9 ]], [], ...]
    base_name_list_1=list()
    base_name_list_2=list()
    base_atom_list_1=list()
    base_atom_list_2=list()

    for i in range(LIST_NUM):

        if os.path.isfile(output_name[i]):
            print "backup %s to %s" %(output_name[i],"#"+output_name[i]+"#")
            try:
                os.rename(output_name[i],"#"+output_name[i]+"#")
            except OSError,e: 
                print e
                print "the file %s will be overwrited!" %output_name[i]

        fp = open(output_name[i], 'w')
        fp.write("#Group 1: ")
        for j in base_list_1[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#Group 2: ")
        for j in base_list_2[i]:
            fp.write("%d\t " %j)
        fp.write("\n")
        fp.write("#skip:%d\n" %skip)
        fp.write("#     L0        L        cos(theta)        theta\n")
        fp.close()

#        base_name_list_1.append( [residue_list[j-1] for j in base_list_1[i]])
#        base_name_list_2.append( [residue_list[j-1] for j in base_list_2[i]])

        temp_list=list()
        for m in base_list_1[i]:
            temp_list.append(residue_list[m-1])
        base_name_list_1.append(temp_list)
        # base_name_list_1.append(residue_list[m-1])

        temp_list=list()
        for m in base_list_2[i]:
            temp_list.append(residue_list[m-1])
        base_name_list_2.append(temp_list)


        base_atom_list_1.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_1[i]])
        base_atom_list_2.append([DNA_matrix.Get_baseID_list(Atom_list,j) for j in base_list_2[i]])
        # print base_atom_list_1
    # print base_atom_list_1, base_atom_list_2

    u=MDAnalysis.Universe(coor_file,traj_file)

    if traj_file.endswith("mdcrd") or traj_file.endswith("dcd"):
        pass
    else:
        try:
            dt=u.trajectory.dt
        except:
            dt=0.0

    for ts in u.trajectory:

        ##construct a list to save the R value to each step
        L_stepList=list()
        R_middle_LIST_1=list()
        R_middle_LIST_2=list()
        O_middle_LIST_1=list()
        O_middle_LIST_2=list()

        time=float((ts.frame-1)*dt)
        print "time %8.1f" %time
        if dt > 0.0:
            if time >= float(begin):
                continue
            if time > float(end) and end !=-1:
                break
        print "time %8.1f" %time
        if ts.frame % skip == 0 :
            print "frame %d" %ts.frame

            #
            #calculate the L value for each base step 
            for i in range(chainLIST_NUM):
                print i+1,
                r1=[]
                '''the group 1 rotate list'''
                r2=[]
                '''the group 2 rotate list'''
                c1=[]
                '''the group 1 coordinate list'''
                c2=[]
                '''the group 2 coordinate list'''
                for m in range(len(chain_bpList_1[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in chain_bpAtomList_1[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), chain_bpList_1[i][m][0])
                    #base_name_list_1[index of the groups][index of the base of group 1][base_name,base_serial]
                    r1.append(result[0])
                    c1.append(result[1])

                for m in range(len(chain_bpList_2[i])):
                    temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in chain_bpAtomList_2[i][m] ]
                    result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), chain_bpList_2[i][m][0])
                    r2.append(result[0])
                    c2.append(result[1])

                middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
                middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])

                R_middle_LIST_1.append(middle_r1)
                R_middle_LIST_2.append(middle_r2)
                O_middle_LIST_1.append(middle_c1)
                O_middle_LIST_2.append(middle_c2)

                a=DNA_param.persistence_length_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
                L_stepList.append(a[1])
            print
            print len(R_middle_LIST_1)==len(R_middle_LIST_2)==len(O_middle_LIST_1)==len(O_middle_LIST_2)==chainLIST_NUM
            print L_stepList

            ####
            
            for i in range(LIST_NUM):
                L0=0.0
                fp = open(output_name[i], 'a')
                print output_name[i]
                bp_INDEX0=0
                bp_INDEX1=0
                print base_name_list_1[0], base_name_list_2[0]
                print chain_bpList_1[1], chain_bpList_2[1]
                for x in range(chainLIST_NUM):
                    if base_name_list_1[i] ==  chain_bpList_1[x]:
                        bp_INDEX0=x
                        continue
                print bp_INDEX0
                for x in range(chainLIST_NUM):
                    if base_name_list_2[i] == chain_bpList_2[x]:
                        bp_INDEX1=x
                        continue
                print bp_INDEX1
                print bp_INDEX0, bp_INDEX1
                middle_r1, middle_c1 = R_middle_LIST_1[bp_INDEX0], O_middle_LIST_1[bp_INDEX0]
                middle_r2,middle_c2 = R_middle_LIST_2[bp_INDEX1], O_middle_LIST_2[bp_INDEX1]

                a=DNA_param.persistence_length_parameters(middle_r1,middle_r2,middle_c1,middle_c2)

                for x in xrange(bp_INDEX0,bp_INDEX1+1):
                    L0+=L_stepList[x]

                fp.write("%8.2f    %8.2f    %8.2f    %8.2f \n" %(L0,a[1],a[2],a[3]))

                if ts.frame % 10 ==0 and i==0:
                    NOW_TIME=Time.time()
                    usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

                fp.close()

            # for i in range(LIST_NUM):
            #     print i
            #     r1=[]
            #     '''the group 1 rotate list'''
            #     r2=[]
            #     '''the group 2 rotate list'''
            #     c1=[]
            #     '''the group 1 coordinate list'''
            #     c2=[]
            #     '''the group 2 coordinate list'''
            #     for m in range(len(base_name_list_1[i])):
            #         temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_1[i][m] ]
            #         result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_1[i][m][0])
            #         #base_name_list_1[index of the groups][index of the base of group 1][base_name,base_serial]
            #         r1.append(result[0])
            #         c1.append(result[1])

            #     for m in range(len(base_name_list_2[i])):
            #         temp_list = [ [ts._x[x-1], ts._y[x-1], ts._z[x-1]] for x in base_atom_list_2[i][m] ]
            #         result = DNA_matrix.Get_rotate_matrix(numpy.array(temp_list), base_name_list_2[i][m][0])
            #         r2.append(result[0])
            #         c2.append(result[1])

            #     fp = open(output_name[i], 'a')
            #     print output_name[i]

            #     # if CALCU=="pair":
            #     #     a=DNA_param.base_pair_parameters(r1[0],r2[0],c1[0],c2[0])
            #     #     fp.write("%8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n" %(a[0],a[1],a[2],a[3],a[4],a[5]))
            #     # else:
            #     middle_r1,middle_c1=DNA_param.middle_frame(r1[0],r1[1],c1[0],c1[1])
            #     middle_r2,middle_c2=DNA_param.middle_frame(r2[0],r2[1],c2[0],c2[1])

            #     a=DNA_param.persistence_length_parameters(middle_r1,middle_r2,middle_c1,middle_c2)
            #     fp.write("%8.2f    %8.2f    %8.2f    %8.2f \n" %(a[0],a[1],a[2],a[3]))

            #     if ts.frame % 10 ==0 and i==0:
            #         NOW_TIME=Time.time()
            #         usage.echo("  analysis frame %6d, time %8.1f ps, time used %8.2f s\r" %(ts.frame, time,NOW_TIME-START_TIME))

            #     fp.close()

    print "The DNA helical analysis finished"
    print "The result are in file: %s" %output_name