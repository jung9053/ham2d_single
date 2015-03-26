 #ifdef UNIT_TEST
 program testcross
 

 for(f=f1;f<f2-idv;f++)
 {
    iface = g->chainConn[f]
    leftCell = g->faces[6*iface+2]  
    rightCell = g->faces[6*iface+4]  
    for(i=0;i<4;i++)
    {
      c1=g->neig[4*leftCell+i]
      if(c1.ne.rightCell)
    }
 }


 end program testcross
 #endif

