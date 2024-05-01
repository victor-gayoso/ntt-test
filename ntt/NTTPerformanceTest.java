package ntt;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author Víctor Gayoso Martínez & David García Lleyda
 */
public class NTTPerformanceTest
{

    long N = 8;
    long Q = 12289;
    long OMEGA = 0;
    long OMEGA_INV = 0;
    long PSI = 0;
    long PSI_INV = 0;
    long N_INV = 0;
    long temp, temp2;
    
    long[][] laPotenciasDoblesOmega;
    long[][] laPotenciasDoblesOmegaInv;  
    long[] laPotenciasOmega;
    long[] laPotenciasOmegaInv;

    long[] laPotenciasPhi;
    long[] laPotenciasPhiInv;  
    
    boolean negacyclic = false;
    boolean debug = true; 
    
    long start, stop, timeElapsedNTTPlain=0, timeElapsedNTTRecursive=0, timeElapsedDirect=0, timeElapsedNTTPlainInit=0, timeElapsedNTTRecursiveInit = 0;
    int iterations = 1;
    
    public static void main(String[] args)
    {
        NTTPerformanceTest  nttpt = new NTTPerformanceTest (args);
    }
    
    public NTTPerformanceTest (String[] args)
    {
        if(args.length == 0)
        {
            System.out.println("INFO:  No input. Default values are used");
            printUsage();
        }
        else
        {

            if(args.length == 3)
            {
                if(args[0].equals("+"))
                {
                    negacyclic = true;
                }
                else
                {
                    if(args[0].equals("-"))
                    {
                        negacyclic = false;
                    }
                    else
                    {
                        System.out.println("\r\nERROR: invalid mode of operation");
                        printUsage();
                        return;                    
                    }
                }
                
                N = Integer.parseInt(args[1]); 
                
                iterations = Integer.parseInt(args[2]); 
                
            }
            else
            {
                printUsage();
                return;
            }
        
        }
    
        ////////////////////////////
        // STEP 1: NTT BASIC INIT
        ////////////////////////////
        
        start = System.nanoTime();
        NTTPlainInit();
        stop = System.nanoTime();
        timeElapsedNTTPlainInit = timeElapsedNTTPlainInit + (stop - start); 
        
        
        ////////////////////////////
        // STEP 2: NTT RECURSIVE INIT        
        ////////////////////////////
        
        start = System.nanoTime();
        NTTRecursiveInit();
        stop = System.nanoTime();
        timeElapsedNTTRecursiveInit = timeElapsedNTTRecursiveInit + (stop - start);         

        if(negacyclic)
        {
            System.out.println("===== PARAMETERS =====\r\n\r\nN: " + N + "\r\nN_INV: "  + N_INV  + "\r\nQ: " + Q + "\r\nOMEGA: " + OMEGA + 
                      "\r\nOMEGA_INV: " + OMEGA_INV + "\r\nPSI: " + PSI + "\r\nPSI_INV: " + PSI_INV  + "\r\nx^(N+1)\r\n");
        }
        else
        {
            System.out.println("===== PARAMETERS =====\r\n\r\nN: " + N + "\r\nN_INV: "  + N_INV + "\r\nQ: " + Q + "\r\nOMEGA: " + OMEGA + 
                      "\r\nOMEGA_INV: " + OMEGA_INV  + "\r\nx^(N-1)\r\n");
        }         
        
        
        //////////////////////////////////////////////
        // STEP 3: GENERATION OF RANDOM COEFFICIENTS
        //////////////////////////////////////////////
        
        long[] a = new long[(int)N];
        long[] b = new long[(int)N];
        long[] aTrans = new long[(int)N]; 
        long[] bTrans = new long[(int)N];
        
        for(int it=0; it < iterations; it++)
        {       
            
            Random random = new Random(System.currentTimeMillis());
        
            for(int i=0; i<N; i++)
            {
                a[i] = random.nextInt((int)Q);
            }
            
            if(debug)
            {
                System.out.println("-----------------------------------------------");               
                System.out.println("Polynomial a(x):                 " + Arrays.toString(a));
            }             
            
            for(int i=0; i<N; i++)
            {
                b[i] = random.nextInt((int)Q);
            }        
            
            if(debug)
            {
                System.out.println("Polynomial b(x):                 " + Arrays.toString(b));         
            }
            
            
            //////////////////////////////////
            // STEP 4: DIRECT METHOD
            //////////////////////////////////
            
            long[] mulDirect = new long[(int)N];
                        
            start = System.nanoTime();
            multiplicationDirect(a,b,mulDirect);
            
            stop = System.nanoTime();
            timeElapsedDirect = timeElapsedDirect + (stop - start);            
            
            
            ////////////////////////////////////
            // STEP 5: NTT BASIC METHOD 
            //////////////////////////////////
            
            start = System.nanoTime();

            NTTPlain(a,aTrans);

            if(debug)
            {
                if(negacyclic)
                {
                    System.out.println("Negacyclic NTT(a(x)):            " + Arrays.toString(aTrans));            
                }
                else
                {
                    System.out.println("Regular NTT(a(x)):               " + Arrays.toString(aTrans));
                }
            }

            NTTPlain(b,bTrans);

            if(debug)
            {
                if(negacyclic)
                {
                    System.out.println("Negacyclic NTT(b(x))):           " + Arrays.toString(bTrans));           
                }
                else
                {
                    System.out.println("Regular NTT(b(x))):              " + Arrays.toString(bTrans));
                }
            }

            long[] mul = new long[(int)N];
            long[] mulTrans = new long[(int)N];

            PointWiseMult(aTrans,bTrans,mulTrans);

            INTTPlain(mulTrans,mul);

            stop = System.nanoTime();
            timeElapsedNTTPlain = timeElapsedNTTPlain + (stop - start);

            for(int i=0; i<N; i++)
            {
                mul[i] = (mul[i]+Q)%Q;
            }
            
            
            ///////////////////////////////////
            // STEP 6: NTT RECURSIVE METHOD
            //////////////////////////////////            
            
            long[] laAmod = new long[(int)N];
            long[] laBmod = new long[(int)N];
            
            start = System.nanoTime();
            
            if(negacyclic)
            {
                for(int i=0; i<N; i++)
                {
                    laAmod[i] = (a[i]*laPotenciasPhi[i])%Q;
                    laBmod[i] = (b[i]*laPotenciasPhi[i])%Q;
                }
            }
            else
            {
                for(int i=0; i<N; i++)
                {                
                    laAmod[i] = a[i];
                    laBmod[i] = b[i];                
                }
            }
            
            aTrans = NTTRecursive(laAmod, laPotenciasOmega);

            bTrans = NTTRecursive(laBmod, laPotenciasOmega);    
            
            PointWiseMult(aTrans,bTrans,mulTrans);
            
            long[] mulNTTrecursive = INTTRecursive(mulTrans,laPotenciasOmegaInv);
            
            for(int i=0; i<N; i++)
            {
                mulNTTrecursive[i] = (mulNTTrecursive[i]*N_INV)%Q;
                
                if(negacyclic)
                {
                   mulNTTrecursive[i] = (mulNTTrecursive[i]*laPotenciasPhiInv[i])%Q;
                }
            }

            stop = System.nanoTime();
            timeElapsedNTTRecursive = timeElapsedNTTRecursive + (stop - start);
            
            if(debug)
            {
                if(negacyclic)
                {
                    System.out.println("Negacyclic NTT(a(x)·b(x)):       " + Arrays.toString(mulTrans));
                } 
                else
                {
                    System.out.println("Regular NTT(a(x)·b(x)):          " + Arrays.toString(mulTrans));
                }
                System.out.println("Direct multiplication a(x)·b(x): " + Arrays.toString(mulDirect)); 
                System.out.println("Plain NTT a(x)·b(x):             " + Arrays.toString(mul)); 
                System.out.println("Recursive NTT a(x)·b(x):         " + Arrays.toString(mulNTTrecursive));
                System.out.println("-----------------------------------------------");  

                for(int i=0; i<(int)N; i++)
                {
                    if(mul[i]!=mulDirect[i])
                    {
                        System.out.println("ERROR: matrices do not match!!!");
                        return;
                    }
                }
            }
        }
        
        System.out.println("END OF COMPUTATION: " + iterations + " iteration/s completed\r\n");

        System.out.println("Timing values provided in nanoseconds\r\n");
        
        System.out.println("-----------------------------------------------");        
        System.out.println("Time NTT Plain Init:              " + (timeElapsedNTTPlainInit));
        System.out.println("Time NTT Recursive Init:          " + (timeElapsedNTTRecursiveInit));        
        System.out.println("Total time Direct:                " + timeElapsedDirect);        
        System.out.println("Total time NTT Plain:             " + timeElapsedNTTPlain);
        System.out.println("Total time NTT Recursive:         " + timeElapsedNTTRecursive);
        System.out.println("-----------------------------------------------");
        System.out.println("Mean time Direct:                 " + (timeElapsedDirect/iterations));
        System.out.println("Mean time NTT Plain:              " + (timeElapsedNTTPlain/iterations));
        System.out.println("Mean time NTT Recursive:          " + (timeElapsedNTTRecursive/iterations));
        System.out.println("-----------------------------------------------");
        System.out.println("Mean time Direct:                 " + (timeElapsedDirect/iterations));
        System.out.println("Mean time NTT Plain + Init:       " + ((timeElapsedNTTPlainInit+timeElapsedNTTPlain)/iterations));
        System.out.println("Mean time NTT Recursive + Init:   " + ((timeElapsedNTTRecursiveInit+timeElapsedNTTRecursive)/iterations));
        System.out.println("-----------------------------------------------");
        System.out.println("Mean time Direct:                 " + (timeElapsedDirect/iterations));
        System.out.println("Mean time 1 NTT Plain + Init:     " + (timeElapsedNTTPlainInit+(timeElapsedNTTPlain/iterations)));
        System.out.println("Mean time 1 NTT Recursive + Init: " + (timeElapsedNTTRecursiveInit+(timeElapsedNTTRecursive/iterations)));
        System.out.println("-----------------------------------------------");
    }
    
    
    void NTTPlainInit()
    {
        BigInteger biN,biQ,biOMEGA,biTemp,biTemp2,biUNO;
        biN = new BigInteger(Integer.toString((int)N));
        biQ = new BigInteger(Integer.toString((int)Q));
        biUNO = new BigInteger("1");
        
        boolean res = false;
        
        ////////////////////////////////////////
        // COMPUTATION OF OMEGA      
        ////////////////////////////////////////
        
        if(OMEGA == 0)
        {
            
            for(int i=2; i< Q; i++)
            {
                biOMEGA = new BigInteger(Integer.toString(i));
                biTemp = biOMEGA.modPow(biN, biQ);                
                
                res = false;

                if(biTemp.compareTo(biUNO) == 0)
                {
                    res = true;

                    for(int j=1; j < N; j++)
                    {
                        biTemp2 = new BigInteger(Integer.toString(j));
                        biTemp = biOMEGA.modPow(biTemp2, biQ);                        
                        
                        if(biTemp.compareTo(biUNO) == 0)
                        {
                            res = false;
                            break;
                        }
                    }

                    if(res)
                    {
                        OMEGA = biOMEGA.intValue();
                        break; 
                    }
                }
            }
            if(!res)
            {
                System.out.println("\r\nERROR: OMEGA NOT FOUND");
                return;
            }
        }

        ////////////////////////////////////////
        // COMPUTATION OF PSI      
        ////////////////////////////////////////

        if(negacyclic)
        {
            for(int i=2; i< Q; i++)
            {
                PSI = i;
                temp = ((long)Math.pow(PSI,2))%Q;

                if(temp == OMEGA)
                {
                    res = true;
                    break;
                }
            }
            
            if(!res)
            {
                System.out.println("\r\nERROR: PSI NOT FOUND");
                return;
            }
        }

        ///////////////////////////////////////////
        // COMPUTATION OF OMEGA_INV FROM OMEGA
        ///////////////////////////////////////////
        
        biTemp = (new BigInteger(Integer.toString((int)OMEGA))).modInverse(new BigInteger(Integer.toString((int)Q)));
        OMEGA_INV = biTemp.intValue();
        
        ////////////////////////////////////////
        // COMPUTATION OF N_INV FROM N      
        ////////////////////////////////////////
        
        biTemp = (new BigInteger(Integer.toString((int)N))).modInverse(new BigInteger(Integer.toString((int)Q)));
        N_INV = biTemp.intValue();   
        
        ////////////////////////////////////////
        // COMPUTATION OF PSI_INV FROM PSI  
        ////////////////////////////////////////
        
        if(negacyclic)
        {
            biTemp = (new BigInteger(Integer.toString((int)PSI))).modInverse(new BigInteger(Integer.toString((int)Q)));
            PSI_INV = biTemp.intValue();           
        }
        
        ////////////////////////////////////////
        // COMPUTATION OF OMEGA DOUBLE MATRIX
        ////////////////////////////////////////
        
        laPotenciasDoblesOmega = new long[(int)N][(int)N];

        for(int j=0; j<N;j++)
        {
            laPotenciasDoblesOmega[0][j]=1;
        }
        
        long temp = 1;
        
        for(int j=0; j<N; j++)
        {
            laPotenciasDoblesOmega[1][j]=temp;
            temp = (temp*OMEGA) % Q;
        }
        
        for(int i=2; i<N; i++)
        {
            for(int j=0; j < N; j++)
            {
                laPotenciasDoblesOmega[i][j]=(laPotenciasDoblesOmega[i-1][j]*laPotenciasDoblesOmega[1][j])%Q;
            }
        }
        
        laPotenciasOmega = new long[(int)N];
        
        for(int i=0; i < N; i++)
        {
            laPotenciasOmega[i]=laPotenciasDoblesOmega[i][1];
        }
        
        ////////////////////////////////////////////
        // COMPUTATION OF OMEGA_INV DOUBLE MATRIX 
        ////////////////////////////////////////////
        
        laPotenciasDoblesOmegaInv = new long[(int)N][(int)N];  
        
        temp = 1;
        
        for(int j=0; j<N;j++)
        {
            laPotenciasDoblesOmegaInv[0][j]=1;
        }
        
        
        for(int j=0; j<N; j++)
        {
            laPotenciasDoblesOmegaInv[1][j]=temp;
            temp = (temp*OMEGA_INV) % Q;
        }   

        for(int i=2; i<N; i++)
        {
            for(int j=0; j < N; j++)
            {
                laPotenciasDoblesOmegaInv[i][j]=(laPotenciasDoblesOmegaInv[i-1][j]*laPotenciasDoblesOmegaInv[1][j])%Q;
            }
        } 
        
        laPotenciasOmegaInv = new long[(int)N];
        
        for(int i=0; i < N; i++)
        {
            laPotenciasOmegaInv[i]=laPotenciasDoblesOmegaInv[i][1];
        }

        ///////////////////////////////////////
        // COMPUTATION OF PSI MATRIX
        ///////////////////////////////////////
        
        if(negacyclic)
        {
            laPotenciasPhi = new long[(int)N];

            laPotenciasPhi[0]=1;
            temp = 1;

            for(int i=1; i<N;i++)
            {
                temp = (temp*PSI)%Q;
                laPotenciasPhi[i]=temp;
            }
        }
        
        //////////////////////////////////////////
        // COMPUTATION OF PSI_INV MATRIX
        //////////////////////////////////////////

        if(negacyclic)
        {
            laPotenciasPhiInv = new long[(int)N];     

            laPotenciasPhiInv[0]=1;
            temp = 1;

            for(int i=1; i<N;i++)
            {
                temp = (temp*PSI_INV)%Q;
                laPotenciasPhiInv[i]=temp;
            }
        }
    }

    
    void NTTRecursiveInit()
    {
        BigInteger biN,biQ,biOMEGA,biTemp,biTemp2,biUNO;
        biN = new BigInteger(Integer.toString((int)N));
        biQ = new BigInteger(Integer.toString((int)Q));
        biUNO = new BigInteger("1");
        
        boolean res = false;
        
        ////////////////////////////////////////
        // COMPUTATION OF OMEGA     
        ////////////////////////////////////////
        
        if(OMEGA == 0)
        {
            for(int i=2; i< Q; i++)
            {
                biOMEGA = new BigInteger(Integer.toString(i));
                biTemp = biOMEGA.modPow(biN, biQ);                
                res = false;

                if(biTemp.compareTo(biUNO) == 0)
                {
                    res = true;

                    for(int j=1; j < N; j++)
                    {
                        biTemp2 = new BigInteger(Integer.toString(j));
                        biTemp = biOMEGA.modPow(biTemp2, biQ);                        
                        
                        if(biTemp.compareTo(biUNO) == 0)
                        {
                            res = false;
                            break;
                        }
                    }

                    if(res)
                    {
                        OMEGA = biOMEGA.intValue();
                        break; 
                    }
                }
            }
            if(!res)
            {
                System.out.println("\r\nERROR: OMEGA NOT FOUND");
                return;
            }
        } 

        ////////////////////////////////////////
        // COMPUTATION OF PSI      
        ////////////////////////////////////////

        if(negacyclic)
        {
            for(int i=2; i< Q; i++)
            {
                PSI = i;
                temp = ((long)Math.pow(PSI,2))%Q;

                if(temp == OMEGA)
                {
                    res = true;
                    break;
                }
            }
            
            if(!res)
            {
                System.out.println("\r\nERROR: PSI NOT FOUND");
                return;
            }
        }

        ///////////////////////////////////////////
        // COMPUTATION OF OMEGA_INV FROM OMEGA
        ///////////////////////////////////////////
        
        biTemp = (new BigInteger(Integer.toString((int)OMEGA))).modInverse(new BigInteger(Integer.toString((int)Q)));
        OMEGA_INV = biTemp.intValue();
        
        ////////////////////////////////////////
        // COMPUTATION OF N_INV FROM N      
        ////////////////////////////////////////
        
        biTemp = (new BigInteger(Integer.toString((int)N))).modInverse(new BigInteger(Integer.toString((int)Q)));
        N_INV = biTemp.intValue();   
        
        ////////////////////////////////////////
        // COMPUTATION OF PSI_INV FROM PSI  
        ////////////////////////////////////////
        
        if(negacyclic)
        {
            biTemp = (new BigInteger(Integer.toString((int)PSI))).modInverse(new BigInteger(Integer.toString((int)Q)));
            PSI_INV = biTemp.intValue();           
        }
        
        ////////////////////////////////////////
        // COMPUTATION OF OMEGA ARRAY
        ////////////////////////////////////////
        
        laPotenciasOmega = new long[(int)N];

        laPotenciasOmega[0]=1;

        for(int i=1; i < N; i++)
        {
            laPotenciasOmega[i]=(laPotenciasOmega[i-1]*OMEGA)%Q;
        }
        
        laPotenciasOmegaInv = new long[(int)N];
        
        laPotenciasOmegaInv[0]=1;
        
        for(int i=1; i < N; i++)
        {
            laPotenciasOmegaInv[i]=(laPotenciasOmegaInv[i-1]*OMEGA_INV)%Q;
        }

        ///////////////////////////////////////
        // COMPUTATION OF PSI ARRAY
        ///////////////////////////////////////
        
        if(negacyclic)
        {
            laPotenciasPhi = new long[(int)N];

            laPotenciasPhi[0]=1;

            for(int i=1; i<N;i++)
            {
                laPotenciasPhi[i]=(laPotenciasPhi[i-1]*PSI)%Q;
            }
        }
        
        //////////////////////////////////////////
        // COMPUTATION OF PSI_INV ARRAY 
        //////////////////////////////////////////

        if(negacyclic)
        {
            laPotenciasPhiInv = new long[(int)N];     

            laPotenciasPhiInv[0]=1;

            for(int i=1; i<N;i++)
            {
                laPotenciasPhiInv[i]=(laPotenciasPhiInv[i-1]*PSI_INV)%Q;
            }
        }
    }
   
    
void multiplicationDirect(long[] laPoly, long[] lbPoly, long[] laRes)
    {
        long[][] laMulti = new long[(int)N][(int)N];
        int desp = 0;
        //int Q2 = Q*Q;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j < N; j++)
            {
                if(negacyclic && (j>(N-1-desp)))
                {
                    laMulti[i][(j+desp)%(int)N] = -lbPoly[i]*laPoly[j]; 
                }
                else
                {
                    laMulti[i][(j+desp)%(int)N] = lbPoly[i]*laPoly[j];                     
                }
            }
            desp++;
        }
        
        for(int j=0; j<N;j++)
        {
            laRes[j]=0;
            for(int i=0; i<N; i++)
            {
                laRes[j] = laRes[j] + laMulti[i][j];
            }
            laRes[j] = ((laRes[j] % Q)+Q)%Q;
        }
    }

   void NTTPlain(long[] laPoly, long[] laPolyModTrans)
    {
        //System.out.println("nplusone: " + nplusone);
        long[] laPolyMod = new long[(int)N];
        long temp = 0;
        
        if(negacyclic)
        {
        for(int i=0; i<N; i++)
            {
                laPolyMod[i] = (laPoly[i]*laPotenciasPhi[i])%Q;
            } 
        }
        
        for (int i=0; i<N; i++)
        {
            temp = 0;
            
            for(int j=0; j < N; j++)
            {
                if(negacyclic) // Para poder hacer mod (x^n+1)
                {
                    temp = (temp + (laPolyMod[j]*(laPotenciasDoblesOmega[i][j]%Q)))%Q;
                }
                else // Para poder hacer mod (x^n-1)
                {
                    temp = (temp + (laPoly[j]*(laPotenciasDoblesOmega[i][j]%Q)))%Q;
                }
            }
            laPolyModTrans[i] = temp;
        }        
    }
    
    long[] NTTRecursive(long[] laPoly, long[] laOmegas)
    {
        if(laPoly.length ==1)
        {
            return laPoly;
        }        
        
        long[] transform = new long[laPoly.length];
        int mitad = laPoly.length/2;
        
        long[] evenPol = new long[laPoly.length/2];
        
        for(int i=0; i<laPoly.length;i=i+2)
        {
            evenPol[i/2] = laPoly[i]; 
        }
        
        long[] oddPol = new long[laPoly.length/2];
        
        for(int i=0; i<laPoly.length;i=i+2)
        {
            oddPol[i/2] = laPoly[i+1]; 
        }      
        
        long[] evenOmegas = new long[laPoly.length/2];
        
        for(int i=0; i<laPoly.length;i=i+2)
        {
            evenOmegas[i/2] = laOmegas[i]; 
        }
        
        long[] even = NTTRecursive(evenPol, evenOmegas);
        long[] odd = NTTRecursive(oddPol, evenOmegas);        

        for(int k=0; k<mitad; k++)
        {
            long w_odd_k = ((laOmegas[k] * odd[k])+Q)%Q;
            long w_even_k = even[k];
        
            transform[k] = ((w_even_k + w_odd_k)+Q) % Q;
            transform[k+mitad] = ((w_even_k - w_odd_k)+Q) % Q;
        }
        
        return transform;

    }
    
    long[] INTTRecursive(long[] laPoly, long[] laOmegasInv)
    {
        if(laPoly.length ==1)
        {
            return laPoly;
        }        
        
        long[] transform = new long[laPoly.length];
        int mitad = laPoly.length/2;
        
        long[] evenPol = new long[laPoly.length/2];
        
        for(int i=0; i<laPoly.length;i=i+2)
        {
            evenPol[i/2] = laPoly[i]; 
        }
        
        long[] oddPol = new long[laPoly.length/2];
        
        for(int i=0; i<laPoly.length;i=i+2)
        {
            oddPol[i/2] = laPoly[i+1]; 
        }      
        
        long[] evenOmegas = new long[laPoly.length/2];
        
        for(int i=0; i<laPoly.length;i=i+2)
        {
            evenOmegas[i/2] = laOmegasInv[i]; 
        }
        
        long[] even = INTTRecursive(evenPol, evenOmegas);
        long[] odd = INTTRecursive(oddPol, evenOmegas);        

        for(int k=0; k<mitad; k++)
        {
            long w_odd_k = ((laOmegasInv[k] * odd[k])+Q)%Q;
            long w_even_k = even[k];
        
            transform[k] = ((w_even_k + w_odd_k)+Q) % Q;
            transform[k+mitad] = ((w_even_k - w_odd_k)+Q) % Q;
        }
        
        return transform;

    }
        
    
    void INTTPlain(long[] laPolyModTrans, long[] laPoly)
    {
        long[] laPolyMod = new long[(int)N];
        for(int i=0; i<N; i++)
        {
            laPolyMod[i]=0;
        }
        long temp = 0;
        
        for (int i=0; i<N; i++)
        {
            temp = 0;
            
            for(int j=0; j < N; j++)
            {
                temp = (temp + (laPolyModTrans[j]*(laPotenciasDoblesOmegaInv[i][j]%Q)))%Q;
            }
            //System.out.println("Suma: " + temp);
            laPolyMod[i] = (temp*N_INV)%Q;
        }   

        for(int i=0; i<N; i++)
        {
            if(negacyclic)
            {
               laPoly[i] = (laPolyMod[i]*laPotenciasPhiInv[i])%Q;
            }
            else
            {
                laPoly[i] = laPolyMod[i];
            }
        }
    }    
              
    void PointWiseMult(long[] laInputA, long[] laInputB, long[] laOutput)
    {
        for(int i=0; i<N; i++)
        {
            laOutput[i] = (laInputA[i]*laInputB[i])%Q;
        }                
    }
    
    void printUsage()
    {
        System.out.println("\r\nUSAGE: java -jar NTTPerformanceTest.jar [+|-] [N] [I]");
        System.out.println("USAGE:   +: negacyclic convolution");
        System.out.println("USAGE:   -: regular convolution");
        System.out.println("USAGE:   N: number of coefficients"); 
        System.out.println("USAGE:   I: number of iterations\r\n");         
    }
}
