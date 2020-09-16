import sys
import math
def extraction_atome(fichier_pdb):
    """ Lit un fichier pdb et renvoie les coordonées de tous les atomes
        impliqués dans la liaison peptidique (C,O,H,N).

    Parametres
    ----------
    fichier_pdb: fichier pdb contenant les coordonnées des atomes
    fichier_text: fichier de sortie qui contient que les coordonnées
                des atomes impliqués dans la liaison hydrogene

    Returns
    -------
    position: liste du numero des differents residus
    fichier_text

    """
    nom = ["C","O","N","H"]
    with open(fichier_pdb, "r") as pdb:
        dico={}
        for ligne in pdb:
            if ligne.startswith("ATOM") and ligne[12:16].strip() in nom:
                posi=ligne[22:27].strip()
                at=ligne[12:16].strip()
                res=ligne[17:20].strip()
                x = ligne[32:38]
                y = ligne[40:46]
                z = ligne[48:54]
                if posi in dico:
                    dico[posi].append((res, at, x, y, z))
                else:
                    dico[posi] = [(res, at, x, y, z)]
    return dico
               

def distance(res1, res2):
    x1, y1, z1=float(res1[0]), float(res1[1]), float(res1[2])
    x2, y2, z2=float(res2[0]), float(res2[1]), float(res2[2])
    norme = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    distance = math.sqrt(norme)
    return distance


def calcul_energie(dico):
    cle=list(dico.keys())
    val=list(dico.values())
    q1 = 0.42
    q2 = 0.28
    f = 332
    liaison = {}
    for i in range(len(cle)):
        res1=dico[cle[i]]
        #print(res1)
        C=res1[1][2:]
        O=res1[2][2:]
        
        for j in range(len(cle)):
            if i!=j:        
                res2 = dico[cle[j]]
                N = res2[0][2:]
            #print((res2))
                if len(res2)==4:
                    H = res2[3][2:]
                    ron=distance(O,N)
                    rch=distance(C,H)
                    roh=distance(O,H)
                    rcn=distance(C,N)
                    energie=round(q1*q2*f*(1/ron + 1/rch - 1/roh - 1/rcn),2)
                    if energie <-0.5:
                        liaison[(cle[i],cle[j])] = energie
    return liaison





if __name__ == '__main__':
    file_pdb = sys.argv[1]
    dico=extraction_atome(file_pdb)
    lh=calcul_energie(dico)
    print(len(lh))


























