import sys
def extraction_atome(fichier_pdb, fichier_text):
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
    modele = 1
    atome = 0
    nom = ["C","O","N","H"]
    with open(fichier_pdb, "r") as pdb, open(fichier_text, "w") as text:
        for ligne in pdb:
            if ligne.startswith("ATOM") and ligne[12:16].strip() in nom:
                atome += 1
                posi=ligne[24:27].strip()
                nam=ligne[12:16].strip()
                resid=ligne[17:20].strip()
                x = ligne[32:38]
                y = ligne[40:46]
                z = ligne[48:54]
                text.write("{}	{}	{}	{}	{}	{} \n".format(posi,resid,nam,x, y, z))
    return fichier_text


def dico_coord(file_text):
	dico = {}
	with open(file_text,"r") as t:
		for ligne in t:
			l = ligne.split("\t")
			posi = l[0]
			res = l[1]
			at = l[2]
			x, y, z = l[3], l[4], l[5].strip()
			if posi in dico:
				dico[posi].append((res, at, x,y,z))
				dico[posi] = [(res, at, x, y, z)]
	return dico, position

if __name__ == '__main__':
    file_pdb = sys.argv[1]
    file_text = sys.argv[2]
    extraction_atome(file_pdb,file_text)
    dico, position=dico_coord(file_text)
    print(position[:10])    

