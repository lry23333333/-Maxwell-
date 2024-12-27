#! /usr/bin/env python 2.7 (62211)
#coding=utf-8

import numpy as np
import time
from abaqus import *
from abaqusConstants import *
from itertools import *

class cohesive_insert():
    skip = 100000000
    nodes_states = {}
    elements_num = 0
    converttomeshpart = 0
    interface_nodes_label_dir = {}
    interface_nodes_label_lis = []

    def __init__(self, Part, Sets):
        #初始参数定义
        self.part = Part
        self.sets = Sets.split(',')
        vpName = session.currentViewportName
        self.model = session.sessionState[vpName]['modelName']
        self.max_nodes_length = len(str(len(mdb.models[self.model].parts[self.part].nodes)))
        self.nodeskip = 10 ** self.max_nodes_length
        

    def convertToMeshpart(self):
    #复制一个初始part，命名为“。。。part-original”，定义插入单元的数量，后续操作将在初始part上进行    
        p = mdb.models[self.model].parts[self.part]
        #初始Part赋给P
        try:
            p.PartFromMesh(name=self.part + '-mesh', copySets=True)
            #复制一个初始part，命名为“初始part名-mesh”
        except:
            print 'error'

        mdb.models[self.model].parts.changeKey(fromName=self.part, toName=self.part + '-original')
        #将初始part更名为“初始part名-original”
        mdb.models[self.model].parts.changeKey(fromName=self.part + '-mesh', toName=self.part)
        #将复制的初始part改回原名
        p = mdb.models[self.model].parts[self.part]
        #重新将初始part赋给p
        cohesive_insert.elements_num = len(p.elements)
        #定义插入单元数量为初始单元数量

    def copyNodes_2D_Multiple(self):
        p = mdb.models[self.model].parts[self.part]
        unionsetName = 'Unionset'
        #节点合集命名
        setNames = self.sets[0].split(';')
        #setNames为初始 人为划分的set名字 的列表
        nodeslabelset_dir = {}
        for i in setNames:
            nodeslabelset_dir[i] = set([ node.label for node in p.sets[i].nodes ])
        #将每个初始set的节点信息提取到字典nodeslabelset_dir{}，key为每个初始set名字    

        sets = [ p.sets[i] for i in setNames ]
        #sets为包括所有 初始set信息 的列表 
        p.SetByBoolean(name=unionsetName, sets=tuple(sets))
        #使用所有初始set信息生成一个联合set
        sets_combinations = list(combinations(setNames, 2))
        #初始set名两两组合

        unionsetName_elements = p.sets[unionsetName].elements
        #存储所有初始单元
        unionsetName_elements_labels = []
        for e in unionsetName_elements:
            unionsetName_elements_labels.append(e.label)
            #存储所有初始单元的标号
        unionsetName_nodes = p.sets[unionsetName].nodes
        #存储所有初始单元的节点
        unionsetName_nodes_labels = []
        for n in unionsetName_nodes:
            unionsetName_nodes_labels.append(n.label)
            #存储所有初始单元节点的标号

        for i in unionsetName_nodes:
            
            # print(i.label)
            i_copy_times = 0

            i_elements = i.getElements()
            #该节点所在的单元
            i_elements_labels = []
            #其所在单元编号列表
            for j in i_elements:
                i_elements_labels.append(j.label)

            for j in i_elements_labels:
                if j in unionsetName_elements_labels:
                    i_copy_times += 1
                    #复制次数为节点所在单元中初始单元的个数

            # for j in i_elements:
            #     if j in unionsetName_elements:
            #         i_copy_times += 1

            if i_copy_times == len(i_elements):
                #如果该节点所在单元全为初始单元
                i_copy_times = i_copy_times - 1
                #复制次数-1
                cohesive_insert.nodes_states[i.label] = {'location': 'edge', 'usecode': 0}
                #若节点所在单元全在原始单元内，则次数-1，节点状态记为edge，使用码记为0
            else:
                cohesive_insert.nodes_states[i.label] = {'location': 'inner', 'usecode': 1}
                #否则记为内部，使用码为1
                #使用过程中这个else并未使用
            for k in range(1, i_copy_times + 1):
                #执行复制次数
                p.Node(i.coordinates, None, k * self.nodeskip + i.label)
                #复制节点，标号为原标号+10**n，坐标不变，次数为其所在初始单元个数-1，即使得最终节点个数为其所在初始单元个数
            print(str(i.label) + '/' + str(len(unionsetName_nodes)))
                

        for combination in sets_combinations:
            #set集的两两组合
            interface_nodes_set = set(nodeslabelset_dir[combination[0]]) & set(nodeslabelset_dir[combination[1]])
            #两初始set节点标号的交集
            if interface_nodes_set != set([]):
                #如果两初始set有接触
                p.SetFromNodeLabels(name='Nset-interface-' + combination[0] + '-' + combination[1], nodeLabels=tuple(interface_nodes_set))
                #通过节点标号的交集，生成两set接触节点的set

        return

    def updateElements_2D_Multiple_TRI3(self):
        p = mdb.models[self.model].parts[self.part]
        cohesiveNodesSetsList = []
        unionsetName = 'Unionset'
        missionElementsLabelSet = []
        nodeslabelset_dir = {}
        elementslabellist_dir = {}
        cohesiveelementslabellists_dir = {}
        setNames = self.sets[0].split(';')
        sets_combinations = list(combinations(setNames, 2))
        interfacenodesunionset = set([])
        interfacenodescombinationlist = []
        for combination in sets_combinations:
            try:
                nodeslabelset_dir['Nset-interface-' + combination[0] + '-' + combination[1]] = set([ node.label for node in p.sets['Nset-interface-' + combination[0] + '-' + combination[1]].nodes ])
                #将两set接触节点的编号存入该字典，key为“'Nset-interface-' + combination[0] + '-' + combination[1]”
                interfacenodescombinationlist.append(combination)
                #将有该初始set组合记入列表
                interfacenodesunionset = interfacenodesunionset | nodeslabelset_dir['Nset-interface-' + combination[0] + '-' + combination[1]]
                #将所有相交节点标号合并，计入数组interfacenodesunionset
            except:
                pass

            cohesiveelementslabellists_dir['Elset-interface-' + combination[0] + '-' + combination[1]] = []

        for setName in setNames:
            nodeslabelset_dir[setName] = set([ node.label for node in p.sets[setName].nodes ])
            #存储初始各个set的节点编号
            elementslabellist_dir[setName] = [ element.label for element in p.sets[setName].elements ]
            #存储初始各个set中的单元编号
            cohesiveelementslabellists_dir[setName] = []


        unionsetName_elements = p.sets[unionsetName].elements
        for i in unionsetName_elements:
        #原始三角单元，点
           # print(i)
            missionElementsLabelSet.append(i.label + cohesive_insert.skip)
            #将初始单元编号提取出来＋skip存入列表
            updatenodes = []
            for j in i.connectivity:
            #j为单元节点列表位置非标号（三角单元由三节点连接而成），非标号
                node_j_label = p.nodes[j].label
                #提取该节点标号到变量
                updatenodes.append(p.nodes.getFromLabel(cohesive_insert.nodes_states[node_j_label]['usecode'] * self.nodeskip + node_j_label))
                #提取该点第'usecode'次复制的节点的标号
                cohesive_insert.nodes_states[node_j_label]['usecode'] += 1
                #为了提取其他初始单元该点的复制节点标号，是一个单元提取出一个节点标号，且与其他单元不同

            nodes = tuple(updatenodes)
            p.Element(nodes=nodes, elemShape=TRI3, label=cohesive_insert.skip + i.label)
            #复制原始三角单元，其节点与其他复制三角单元无交集，标号增加了skip

        p.SetFromElementLabels(name='Elset-' + unionsetName, elementLabels=tuple(missionElementsLabelSet))
        #将所有初始单元复制体存入set——"Elset-'UnionsetNmae'"中，所有复制体无相交节点
        for i in p.sets['Elset-' + unionsetName].elements:
        #复制单元set集
           # print(i)
            i_label = i.label
            print(str(i_label-cohesive_insert.skip) + '/' + str(len(p.sets['Elset-' + unionsetName].elements)))
            i_elem_edges = i.getElemEdges()
            #获得单元的边线
            i_elem_edges_0_nodes = i_elem_edges[0].getNodes()
            #对应边线所在节点位置
            i_elem_edges_1_nodes = i_elem_edges[1].getNodes()
            i_elem_edges_2_nodes = i_elem_edges[2].getNodes()
            m = [
             (
              i_elem_edges_0_nodes[0].label, i_elem_edges_0_nodes[1].label),
             (
              i_elem_edges_1_nodes[0].label, i_elem_edges_1_nodes[1].label),
             (
              i_elem_edges_2_nodes[0].label, i_elem_edges_2_nodes[1].label)]
             #将复制体三角单元边线以节点标号的形式存入m中

            for j in p.elements.getFromLabel(i.label - cohesive_insert.skip).getAdjacentElements():
            #该复制单元所对应初始单元的相邻单元（可确保相邻单元都是初始单元，而不是插入单元）
            #p.sets[Elset-Unionset].elements[0].label
                if j in unionsetName_elements:
                    j_label = j.label
                    #如果该复制体的本体的相邻单元为初始单元
                    j_elem_edges = p.elements.getFromLabel(j.label + cohesive_insert.skip).getElemEdges()
                    #该初始单元的复制体边线
                    j_elem_edges_0_nodes = j_elem_edges[0].getNodes()
                    j_elem_edges_1_nodes = j_elem_edges[1].getNodes()
                    j_elem_edges_2_nodes = j_elem_edges[2].getNodes()
                    n = [
                     (
                      j_elem_edges_0_nodes[0].label,
                      j_elem_edges_0_nodes[1].label),
                     (
                      j_elem_edges_1_nodes[0].label,
                      j_elem_edges_1_nodes[1].label),
                     (
                      j_elem_edges_2_nodes[0].label,
                      j_elem_edges_2_nodes[1].label)]
                     #将复制体相邻的复制体单元边线以节点标号的形式存入n中
                else:
                    #否则就是相邻的复制体单元
                    j_label = j.label
                    j_elem_edges = j.getElemEdges()
                    j_elem_edges_0_nodes = j_elem_edges[0].getNodes()
                    j_elem_edges_1_nodes = j_elem_edges[1].getNodes()
                    j_elem_edges_2_nodes = j_elem_edges[2].getNodes()
                    n = [
                     (
                      j_elem_edges_0_nodes[0].label, j_elem_edges_0_nodes[1].label),
                     (
                      j_elem_edges_1_nodes[0].label, j_elem_edges_1_nodes[1].label),
                     (
                      j_elem_edges_2_nodes[0].label, j_elem_edges_2_nodes[1].label)]
                for k in m:
                    #复制体三角单元的三条边线，节点标号形式
                    for l in n:
                        #复制体相邻单元的三条边线，节点标号形式
                        set_k = set(np.array(k) % self.nodeskip)
                        #求余数，即回归原始节点标号
                        set_l = set(np.array(l) % self.nodeskip)
                        if set_k == set_l and k != l:
                            #如果投影点重合，自身不重合
                            if set_k not in cohesiveNodesSetsList:
                                #确保插入过的线不会再执行插入
                                cohesive_insert.elements_num += 1
                                p.Element(nodes=(
                                 p.nodes.getFromLabel(k[1]), p.nodes.getFromLabel(k[0]), p.nodes.getFromLabel(l[1]),
                                 p.nodes.getFromLabel(l[0])), elemShape=QUAD4, label=cohesive_insert.elements_num)
                                #生成cohesive单元，为两投影重合线连接
                                cohesiveNodesSetsList.append(set_k)
                                ij = (i_label,j_label)
                                set_ij = set(np.array(ij) % cohesive_insert.skip)
                                # print ('set_ij',set_ij)
                                grain = []
                                
                                for setName in setNames:
                                    # print setName
                                    if len(set_ij&set(np.array(elementslabellist_dir[setName]))) == 2:
                                        cohesiveelementslabellists_dir[setName].append(cohesive_insert.elements_num)
                                        break
                                    elif len(set_ij&set(np.array(elementslabellist_dir[setName]))) == 1:
                                        grain.append(setName)
                                        # print('1111')
                                if len(grain)==2:
                                    # print('2222')
                                    for combination in sets_combinations:
                                        if grain[0] in combination and grain[1] in combination:
                                            cohesiveelementslabellists_dir['Elset-interface-' + combination[0] + '-' + combination[1]].append(cohesive_insert.elements_num)
                                            break
                                    break
                                elif len(grain)==1:
                                    print ('error',grain)
                                    break
                                

                                # for combination in sets_combinations:
                                #     if grain[0] in combination and grain[1] in combination:
                                #         cohesiveelementslabellists_dir['Elset-interface-' + combination[0] + '-' + combination[1]].append(cohesive_insert.elements_num)
                                #         grain = []
                                #         break
                                # for setName in setNames:
                                #     for combination in interfacenodescombinationlist:
                                #         d = 0
                                #         if set_k.issubset(nodeslabelset_dir['Nset-interface-' + combination[0] + '-' + combination[1]]):
                                #             d += 1
                                #             break
                                #     #切分集合
                                #     if set_ij.issubset(elementslabellist_dir[setName]):
                                #         cohesiveelementslabellists_dir[setName].append(cohesive_insert.elements_num)
                                #         break
                                #     elif d==1:
                                #         cohesiveelementslabellists_dir['Elset-interface-' + combination[0] + '-' + combination[1]].append(cohesive_insert.elements_num)
                                #         break

                                    #if set_k.issubset(nodeslabelset_dir[setName]) and d == 0:
                                        #A。issubset（B），判断A是否为B的子集
                                    #    cohesiveelementslabellists_dir[setName].append(cohesive_insert.elements_num)
                                    #    break
                                    #elif d==1 and not set_ea.issubset(elementslabellist_dir[setName]):
                                        #cohesiveelementslabellists_dir['Elset-interface-' + combination[0] + '-' + combination[1]].append(cohesive_insert.elements_num)
                                        #break
                                    #else:
                                        #cohesiveelementslabellists_dir[setName].append(cohesive_insert.elements_num)
                                       # break


        num1=0
        num2=0
        num3=0
        num4=0
        for i in missionElementsLabelSet:
            p.deleteElement(p.elements.getFromLabel(i - cohesive_insert.skip))
            #删除该复制体的本体单元
            num1+=1
            print(str(num1)+'/'+str(len(missionElementsLabelSet)))

        for i in p.sets['Elset-' + unionsetName].elements:
            i.setValues(i.label - cohesive_insert.skip)
            #将该复制体单元标号改为本体标号
            num2+=1
            print(str(num2)+'/'+str(len(p.sets['Elset-' + unionsetName].elements)))

        for setName in setNames:
            # print(cohesiveelementslabellists_dir[setName])
            p.SetFromElementLabels(name='Elset-' + setName, elementLabels=tuple(elementslabellist_dir[setName]))
            #建立复制体三角单元set集，标号为原单元标号
            p.SetFromElementLabels(name='Elset-' + setName + '-cohesive', elementLabels=tuple(cohesiveelementslabellists_dir[setName]))
            #建立晶内cohesive单元set集
            num3+=1
            print(str(num3)+'/'+str(len(setNames)))
        for combination in interfacenodescombinationlist:
            p.SetFromElementLabels(name='Elset-interface-' + combination[0] + '-' + combination[1], elementLabels=tuple(cohesiveelementslabellists_dir['Elset-interface-' + combination[0] + '-' + combination[1]]))
            #建立晶面cohesive单元集
            num4+=1
            print(str(num4)+'/'+str(len(interfacenodescombinationlist)))

def mission_execute(Dimension, Part, Mode, Sets, ElementShape):
    t0 = time.time()
    x = cohesive_insert(Part, Sets)
    updateElements_method_dir = {
       (2, 30, 23): x.updateElements_2D_Multiple_TRI3,
       }
       #以（2，30，23）指代函数x.updateElements_2D_Multiple_TRI3
    copyNodes_method_dir = {
       (2, 30): x.copyNodes_2D_Multiple, 
       }
    if cohesive_insert.converttomeshpart == 0:
        x.convertToMeshpart()
        cohesive_insert.converttomeshpart = 1
    # print(Dimension, Mode, ElementShape)
    copyNodes_method_dir[(Dimension, Mode)]()
    #执行cohesive_insert.updateElements_2D_Multiple_TRI3函数
    updateElements_method_dir[(Dimension, Mode, ElementShape)]()
    #执行cohesive_insert.convertToMeshpart()函数
    t1 = time.time()
    print(t1-t0)
    #输出执行时间