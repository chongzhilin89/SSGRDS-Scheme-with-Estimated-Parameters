Data SSGR_DSe_steady;

BigM=50;

mu=0;
sigma=1;
delta=0.2;
n1=5;
n2=10;
L=2.297;
L1=0.991;
L2=1.428;
L3=14;
n=8;
m=30;
array X{*} X1-X100;
array CRLUU{*} CRLUU1-CRLUU10000;
array CRLLL{*} CRLLL1-CRLLL10000;

Do a=1 to 20000;
   Xsumsum=0;Ssquare=0;
   Do k=1 to m;
      Xsum=0;Sumsq=0;
      Do h=1 to n;
         X(h)=mu+sigma*rannor(55555);
		 Xsum=Xsum+X(h);
	  End;
	  Xbar=Xsum/n;
	  Xsumsum=Xsumsum+Xsum;
      Do h=1 to n;
         Sumsq=Sumsq+(X(h)-Xbar)**2;
      End;
      Ssquare=Ssquare+Sumsq; 
   End; 
   sigma0=sqrt(Ssquare/(m*(n-1)));
   mu0=Xsumsum/(m*n);
      
   ok=1;
   i=0;
   CRLU=0;CRLL=0;p=0;obs=0;
   Do while (ok=1);
      i=i+1;

	If (i<BigM+1) then 
    do;
	  Ysum=0;
	  Do j=1 to n1+n2;
         Y=0+1*rannor(33333); /*delta=0, sigma=1 for in-control*/
         Ysum=Ysum+Y;		 
		 /*remove obs=obs+1*/
		 If(j=n1) then do;
		 				  Ybar=Ysum/n1;
						  Z1=(Ybar-mu0)/(sigma0/sqrt(n1));
						  If(-L1<Z1<L1) then do;
												j=700;k=1;CRLU=CRLU+1;CRLL=CRLL+1; /*k=1 marks conforming sample*/
											 end;
						  If(Z1<-L) or (Z1>L) then do;
						                              j=700;p=p+1; 						  							  
	                                                  If(Z1>L) then do;
                                                                       If(k=1) then do;
                                                                                       CRLU=0;  
	                                                                                   CRLUU(p)=CRLU;CRLLL(p)=L3+1; /*CRLUU(p) denotes p_th upper-sided nonconforming CRL value*/ 
																					end;
																			   else do;
												                                       CRLUU(p)=CRLU;CRLLL(p)=L3+1;
													                                end;
						                                            end;
                                                      If(Z1<-L) then do;
                                                                        If(k=1) then do;
                                                                                        CRLL=0;  
	                                                                                    CRLLL(p)=CRLL;CRLUU(p)=L3+1; /*CRLLL(p) denotes p_th lower-sided nonconforming CRL value*/ 
																					 end;
																				else do;
												                                        CRLLL(p)=CRLL;CRLUU(p)=L3+1;
													                                 end;
						                                             end;
													  k=0; /*k=0 marks nonconforming sample*/
									               end;	
					   end;
		 If(j=n1+n2) then do;
		 					 Ybar=Ysum/(n1+n2); 
							 Z=(Ybar-mu0)/(sigma0/sqrt(n1+n2));	
							 If(-L2<Z<L2) then do;
												  k=1;CRLU=CRLU+1;CRLL=CRLL+1;
											   end;
										  else do;										  		  
												  p=p+1;
	                                              If(Z>L2) then do;
                                                                   If(k=1) then do;
                                                                                   CRLU=0;  
	                                                                               CRLUU(p)=CRLU;CRLLL(p)=L3+1; 
																				end;
																		   else do;
																		           CRLUU(p)=CRLU;CRLLL(p)=L3+1;
																				end;
						                                        end;
                                                  If(Z<-L2) then do;
                                                                    If(k=1) then do;
                                                                                    CRLL=0;  
	                                                                                CRLLL(p)=CRLL;CRLUU(p)=L3+1;
																				 end;
																			else do;
																			        CRLLL(p)=CRLL;CRLUU(p)=L3+1;
													                             end;
						                                         end;
												  k=0;
											   end;
						  end;	
	  End;
   end;

   else
   do;
      Ysum=0;
	  Do j=1 to n1+n2;
         Y=delta+sigma*rannor(33333);
         Ysum=Ysum+Y;		 
		 obs=obs+1;
		 If(j=n1) then do;
		 				  Ybar=Ysum/n1;
						  Z1=(Ybar-mu0)/(sigma0/sqrt(n1));
						  If(-L1<Z1<L1) then do;
												j=700;k=1;CRLU=CRLU+1;CRLL=CRLL+1; /*k=1 marks conforming sample*/
											 end;
						  If(Z1<-L) or (Z1>L) then do;
						                              j=700;p=p+1; 						  							  
	                                                  If(Z1>L) then do;
                                                                       If(k=1) then do;
                                                                                       CRLU=CRLU+1;  
	                                                                                   CRLUU(p)=CRLU;CRLLL(p)=L3+1; /*CRLUU(p) denotes p_th upper-sided nonconforming CRL value*/ 
																					end;
																			   else do;
												                                       CRLUU(p)=CRLU;CRLLL(p)=L3+1;
													                                end;
						                                            end;
                                                      If(Z1<-L) then do;
                                                                        If(k=1) then do;
                                                                                        CRLL=CRLL+1;  
	                                                                                    CRLLL(p)=CRLL;CRLUU(p)=L3+1; /*CRLLL(p) denotes p_th lower-sided nonconforming CRL value*/ 
																					 end;
																				else do;
												                                        CRLLL(p)=CRLL;CRLUU(p)=L3+1;
													                                 end;
						                                             end;
													  k=0; /*k=0 marks nonconforming sample*/
									               end;	
					   end;
		 If(j=n1+n2) then do;
		 					 Ybar=Ysum/(n1+n2); 
							 Z=(Ybar-mu0)/(sigma0/sqrt(n1+n2));	
							 If(-L2<Z<L2) then do;
												  k=1;CRLU=CRLU+1;CRLL=CRLL+1;
											   end;
										  else do;										  		  
												  p=p+1;
	                                              If(Z>L2) then do;
                                                                   If(k=1) then do;
                                                                                   CRLU=CRLU+1;  
	                                                                               CRLUU(p)=CRLU;CRLLL(p)=L3+1; 
																				end;
																		   else do;
																		           CRLUU(p)=CRLU;CRLLL(p)=L3+1;
																				end;
						                                        end;
                                                  If(Z<-L2) then do;
                                                                    If(k=1) then do;
                                                                                    CRLL=CRLL+1;  
	                                                                                CRLLL(p)=CRLL;CRLUU(p)=L3+1;
																				 end;
																			else do;
																			        CRLLL(p)=CRLL;CRLUU(p)=L3+1;
													                             end;
						                                         end;
												  k=0;
											   end;
						  end;	
	  End;
	  If(k=0) then do;
                      If(p=1 and CRLUU(p)<=L3) then do;
                                                       ARL=i-BigM;ANOS=obs;
                                                       output;ok=0;
									                end;

                      If(p=1 and CRLLL(p)<=L3) then do;
                                                       ARL=i-BigM;ANOS=obs;
                                                       output;ok=0;
									                end;

                      If(p>=3 and CRLUU(p-1)<=L3 and CRLUU(p)<=L3) then do;
                                                                           ARL=i-BigM;ANOS=obs;
                                                                           output;ok=0;
									                                    end;

                      If(p>=3 and CRLLL(p-1)<=L3 and CRLLL(p)<=L3) then do;
                                                                           ARL=i-BigM;ANOS=obs;
                                                                           output;ok=0;
									                                    end;
                      CRLU=0;CRLL=0;
				   end;
	end;
   End;
End;

run;
proc means;
var ARL ANOS;
run;
