function Basins_Thule=GenMask_Basins_Thule_1km(Basins_Thule)
    %Heilprin Gl.
     Basins_Thule(Basins_Thule==233|Basins_Thule==878|Basins_Thule==1336|Basins_Thule==1019)=-12;
    temp=Basins_Thule(187:213, 134:147);
    temp(temp~=-12) = -12;
	Basins_Thule(187:213, 134:147)=temp;
    %Tracy Gl.
    Basins_Thule(Basins_Thule==286|Basins_Thule==118|Basins_Thule==1399)=-13;
    %small one 
    Basins_Thule(Basins_Thule==405)=-14;
    Basins_Thule(Basins_Thule==438)=-15;
    Basins_Thule(Basins_Thule==710)=-16;
    
    Basins_Thule(Basins_Thule==258|Basins_Thule==451)=-17;
    Basins_Thule(Basins_Thule==293)=-18;
    Basins_Thule(Basins_Thule==114)=-19;
    Basins_Thule(Basins_Thule==351)=-20;
    Basins_Thule(Basins_Thule==1027|Basins_Thule==1764)=-21;
    Basins_Thule(Basins_Thule==177)=-22;
    Basins_Thule(Basins_Thule==25|Basins_Thule==119|Basins_Thule==55)=-23;
end