cutoff=3.0;
% load protein and ligand pdb files
protein=pdbread('protein.pdb');
ligand=pdbread('ligand.pdb');
% calculating protein coordinates,and adding the radius of each atom present in protein
[s1,s2]=size(protein.Model.Atom);
for n=1:s2
    proteinXYZ(n,1:3)=[protein.Model.Atom(n).X protein.Model.Atom(n).Y protein.Model.Atom(n).Z];
    proteinAtom=protein.Model.Atom(n).element;
    if proteinAtom=='H'
        r=0.25;
    elseif proteinAtom=='C'
        r=0.70;
    elseif proteinAtom=='N'
        r=0.65;
    elseif proteinAtom=='O'
        r=0.60;
    elseif proteinAtom=='S'
        r=1.00;
    elseif proteinAtom=='F'
        r=0.50;
    elseif proteinAtom=='CL'
        r=1.00;
    elseif proteinAtom=='P'
        r=1.00;
    elseif proteinAtom=='BR'
        r=1.15;
    elseif proteinAtom=='I'
        r=1.40;
    end
    proteinR(n,1)=r;
end
% calculating ligand coordinates,and adding the radius of each atom present in ligand
[s3,s4]=size(ligand.Model.Atom);
for n=1:s4;
    ligandXYZ(n,1:3)=[ligand.Model.Atom(n).X ligand.Model.Atom(n).Y ligand.Model.Atom(n).Z];
    ligandAtom=ligand.Model.Atom(n).element;
    if ligandAtom=='H'
        r=0.25;
    elseif ligandAtom=='C'
        r=0.70;
    elseif ligandAtom=='N'
        r=0.65;
    elseif ligandAtom=='O'
        r=0.60;
    elseif ligandAtom=='S'
        r=1.00;
    elseif ligandAtom=='F'
        r=0.50;
    elseif ligandAtom=='CL'
        r=1.00;
    elseif ligandAtom=='P'
        r=1.00;
    elseif ligandAtom=='BR'
        r=1.15;
    elseif ligandAtom=='I'
        r=1.40;
    end
    ligandR(n,1)=r;
end
% calculating distances between protein-ligand atom pairs
for n=1:s4;
    for nn=1:s2;
        d(n,nn)=sqrt(sum((ligandXYZ(n,:)-proteinXYZ(nn,:)).^2))-ligandR(n,1)-proteinR(nn,1);
    end
end
% calculating DB_Score with specific cutoff
for n=1:s4
    for nn=1:s2
        if d(n,nn)>cutoff
            d(n,nn)=0;
        end
    end
end
q=(max(d'))';
s=0;
for n=1:s4
    if q(n,1)~=0
        s=s+sum(d(n,:));
    end
end
DB_score=s
