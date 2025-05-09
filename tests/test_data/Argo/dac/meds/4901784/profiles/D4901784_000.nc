CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  
   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       MEDS   source        
Argo float     history       2018-12-20T17:00:49Z creation      
references        (http://www.argodatamgt.org/Documentation   comment              user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         IPRIMARY|https://orcid.org/0000-0002-1716-6352|Zhimin(Robert) Ma, OSB, DFO         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8$   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    84   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8D   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8L   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9    DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9(   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9,   PLATFORM_TYPE                     	long_name         Type of float      
_FillValue               conventions       Argo reference table 23          90   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9P   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9p   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T      
resolution        >�����h�        9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�����h�        9�   LATITUDE            	   	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y      	reference         WGS84      coordinate_reference_frame        urn:ogc:crs:EPSG::4326          9�   	LONGITUDE               	   	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X      	reference         WGS84      coordinate_reference_frame        urn:ogc:crs:EPSG::4326          9�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :�   PRES         
         	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z      coordinate_reference_frame        urn:ogc:crs:EPSG::5113       (  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   standard_name         sea_water_pressure     axis      X        (  E   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   M0   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     (  O<   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     (  Wd   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   _�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o   standard_name         sea_water_temperature        (  a�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   i�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     (  k�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     (  s�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   |   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o   standard_name         sea_water_salinity       (  ~(   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   �P   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     (  �\   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    
_FillValue               conventions       YYYYMMDDHHMISS        T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �8   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �<   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �@   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �D   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �H   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��Argo profile    3.1 1.2 19500101000000  20181220120050  20210315164152  4901784 Argo Canada                                                     Blair Greenan                                                   PRES            TEMP            PSAL                A   ME  4901784_9979_PF                 2C+ D   NOVA                            200                             n/a                             865 @�]�""""1   @�]�""""@G?��   �` ;�   1   GPS     A   A   A   Primary sampling: discrete                                                                                                                                                                                                                                         @L��@�33@�  @���A33A.ffAI��AfffA���A�ffA�  A���A���A�  A���A�33A�  A�  A�  BffB
ffB��B��BffB&��B-��B2��B8ffB>ffBDffBJ��BP��BW33B]��Bd  BjffBq33Bw��B~ffB�  B���B�33B���B�ffB�33B�  B���B���B���B�ffB�33B�  B���B�  B���B�33B���B�ffB�  B̙�B�33B�33B�ffB�ffB�33B�ffB�33B�B�33B���C�3C��C33C	�3C33C�3C33C  C��CL�C��C�C �C"ffC$��C'�fC*L�C,�3C/  C1� C4  C6ffC8�fC;� C>  C@ffCC  CE� CH  CJ� CM�CO��CR�CT�3CW33CY��C\ffC^�C`�3CcffCe�fCh��CkL�Cn  Co��Cr� CuL�Cw33Cz  C|��C~�3C�� C�&fC��C���C��3C��fC�Y�C�� C��3C��3C��C���C�� C��3C�Y�C�L�C��3C��C�  C�ffC���C�� C�&fC���C�� C��3C�Y�C�Y�C�� C�&fC�&fC��C��fC��C�Y�C���C�� C�  C�33C�ffC��3C��fC�33C�s3C��fC��fC�&fC�ffC��fC��fC�&fC�s3C��3C�  C�@ CŌ�C���C��C�Y�CʦfCˀ C�ٚC��C�s3C���C��C�ffC�L�CզfC�  C�Y�Cٳ3Cڙ�C��3C�L�C�@ Cߙ�C��3C�Y�C�3C䙚C�  C�ffC�Y�C�� C�&fC��C� C��3C��fC�Y�C���C�� C�� C�&fC��fC�ٚC��C�� C�Y�D �fD  D` D�fD�fD&fD` D	� D
ٚD�D` D� D� D&fDl�D�3D�fDٚD&fDs3D� D�D33DS3D�3D 3D!9�D"y�D#�fD$��D&&fD'Y�D(�3D*  D+l�D,��D-�fD/&fD0ffD1�fD2�fD4,�D5@ D6S3D7� D8��D:@ D;�3D<� D>,�D?FfD@ffDA��DC3DD33DES3DF�3DH�DI,�DJY�DK��DM  DNFfDOs3DP� DQ��DS&fDT` DU��DV�3DX3DYS3DZ�3D[�3D]3D^Y�D_�fD`�3Db9�Dc� Dd�3Df&fDg9�DhS3Di�fDj��DlS3Dms3Dn�3Do��DqL�Dr��Ds��Dt�3DvS3Dw�fDxٚDz3D{L�D|��D}��D�D�0 D��3D�\�D��D���D�33D�ٚD�� D�)�D��3D�\�D���D���D�C3D��3D��3D�fD��fD�Y�D��fD���D�6fD��3D�p D� D�� D�P D�� D�� D�0 D���D�� D�  D��fD�ffD���D��3D�FfD���D�vfD�,�D��fD�\�D��fD�� D�FfD���D�vfD�0 D�ɚD�c3D���D���D�6fD�� D���D�)�D�ɚD�ffD��3D��3D�P D�� D�p D�  D�� D�` D�	�D���D�C3D�� D�|�D��D���D�\�D�  D�� D�@ D�� D��3D�#3D��fD�i�D� D��fD�@ D���D�s3D��D��fD�l�D���D���D�9�D��D�y�D�#3D��fD�ffD�3D��3D�@ D�� D��3D�&fD���D�s3D�  DÉ�D�0 D�ٚDŉ�D��DƬ�D�\�D�  DȠ D�VfD��3Dʐ D�,�D���D�l�D��DͰ D�S3D��fDπ D�fDЬ�D�S3D���DҦfD�P D�� D�p D��D���D�` D�� Dנ D�P D��3D�y�D�)�D�� D�S3D�  Dܬ�D�I�D��Dމ�D�)�D�ɚD�i�D��D��D�P D��3D�fD�  D��D�P D��fD� D�I�D��3D�|�D��D�fD�c3D�3D�3D�6fD��fD�fD�,�D��D�Y�D��fD�fD�33D��3D�s3D�3D�3D�S3D��fD��fD�6fD�ٚD�y�D��D���D�ffD�3D��fD�VfD�  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @S33@�ff@�33@�  A��A0  AK34Ah  A���A�33A���A�fgA���A���A͙�A�  A���A���A���B��B
��B33B  B��B'33B.  B333B8��B>��BD��BK33BQ33BW��B^  BdffBj��Bq��Bx  B~��B�33B���B�ffB�  B���B�ffB�33B�  B���B���B���B�ffB�33B���B�33B���B�ffB�  Bř�B�33B���B�ffB�ffBۙ�BᙙB�ffB뙙B�ffB���B�ffB���C��C�gCL�C	��CL�C��CL�C�C�4CfgC�4C34C 34C"� C$�gC(  C*fgC,��C/�C1��C4�C6� C9  C;��C>�C@� CC�CE��CH�CJ��CM34CO�4CR34CT��CWL�CY�gC\� C^34C`��Cc� Cf  Ch�4CkfgCn�Co�gCr��CufgCwL�Cz�C|�gC~��C���C�33C�&gC���C�  C��3C�fgC���C�� C�� C�&gC���C���C�  C�fgC�Y�C�� C�&gC��C�s3C�ٚC���C�33C���C���C�  C�fgC�fgC���C�33C�33C�&gC��3C�&gC�fgC���C���C��C�@ C�s3C�� C��3C�@ C�� C��3C��3C�33C�s3C��3C��3C�33C�� C�� C��C�L�Cř�C�ٚC�&gC�fgCʳ3Cˌ�C��gC�&gCπ C�ٚC�&gC�s3C�Y�Cճ3C��C�fgC�� CڦgC�  C�Y�C�L�CߦgC�  C�fgC�� C�gC��C�s3C�fgC���C�33C�&gC��C�  C��3C�fgC�ٚC���C���C�33C��3C��gC��C���C�fgD ��D&fDffD��D��D,�DffD	�fD
� D  DffD�fD�fD,�Ds3D��D��D� D,�Dy�D�fD  D9�DY�D��D �D!@ D"� D#��D%  D&,�D'` D(��D*fD+s3D,�3D-��D/,�D0l�D1��D2��D433D5FfD6Y�D7�fD8�3D:FfD;��D<�fD>33D?L�D@l�DA� DC�DD9�DEY�DF��DH3DI33DJ` DK� DM&fDNL�DOy�DP�fDQ�3DS,�DTffDU� DVٙDX�DYY�DZ��D[ٙD]�D^` D_��D`��Db@ Dc�fDdٙDf,�Dg@ DhY�Di��Dk  DlY�Dmy�Dn��Do�3DqS3Dr�3Ds�3Dt��DvY�Dw��Dx� Dz�D{S3D|�3D}�3D  D�33D��fD�` D���D�� D�6fD���D��3D�,�D��fD�` D�� D���D�FfD��fD��fD��D�əD�\�D���D�� D�9�D��fD�s3D�3D��3D�S3D��3D��3D�33D�� D��3D�#3D���D�i�D�  D��fD�I�D�� D�y�D�0 D�əD�` D���D��3D�I�D�� D�y�D�33D���D�ffD�  D���D�9�D��3D�� D�,�D���D�i�D��fD��fD�S3D��3D�s3D�#3D��3D�c3D��D���D�FfD��3D�� D��D���D�` D�3D��3D�C3D��3D��fD�&fD�əD�l�D�3D���D�C3D�� D�vfD�  D�əD�p D�  D�� D�<�D���D�|�D�&fD�əD�i�D�fD��fD�C3D��3D��fD�)�D�� D�vfD�3DÌ�D�33D���DŌ�D��Dư D�` D�3Dȣ3D�Y�D��fDʓ3D�0 D�� D�p D� Dͳ3D�VfD���Dσ3D�	�Dа D�VfD�  Dҩ�D�S3D��3D�s3D�  D�� D�c3D��3Dף3D�S3D��fD�|�D�,�D��3D�VfD�3Dܰ D�L�D���Dތ�D�,�D���D�l�D� D� D�S3D��fD㙙D�#3D� D�S3D���D�3D�L�D��fD� D� D鹙D�ffD�fD�fD�9�D��D홙D�0 D�� D�\�D���D�D�6fD��fD�vfD�fD�fD�VfD���D���D�9�D���D�|�D�  D�� D�i�D�fD���D�Y�D�3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�v�A�~�A�~�A�z�A�jA�p�A�~�A�hsA��A�p�A�bNA�I�A�I�A�{A~r�Av5?Aex�Ac�-Ab��A^^5AW+ATv�AQ"�AO/AL-AJ$�AI\)AH�AEAD�ABVAA��A>^5A=�;A=�#A=��A<�jA;`BA:~�A9�A8��A8�DA7�wA7"�A77LA5��A5�A61'A6(�A5C�A5A4�A3G�A/��A-��A,�jA+��A+S�A(JA$~�A"VA�^A33A��A�A��A�FAhsA�yAC�@�\)@��T@�z�@��@�{@�?}@��@�j@�x�@��@�`B@�/@�bN@��@�\)@�@�@��@�bN@���@��m@�\)@�@��@���@�Ĝ@�Q�@�@�A�@��H@�h@�r�@�ȴ@��@�(�@�
=@�@�V@��;@�\@�@�dZ@�M�@��T@�%@�ƨ@�"�@�E�@�1'@ա�@��/@�Q�@ӕ�@���@ѡ�@� �@Ώ\@�X@˾w@�K�@�M�@Ȭ@�C�@Ɨ�@��`@�bN@���@���@�I�@��@���@���@���@��@�n�@�5?@�p�@��@�dZ@��y@�$�@�X@�Ĝ@���@�5?@��@�{@��-@��;@�"�@�^5@��@�;d@�@�(�@��#@�n�@��T@��7@�9X@�@��@�=q@�
=@��@��y@�A�@��@��@�5?@�ff@��@��@��@���@�j@���@��@���@�33@��j@���@���@��#@�b@�O�@�p�@�V@��@�  @��@�Ĝ@��@��@��#@�M�@�$�@�`B@�x�@�%@�r�@���@���@��m@���@���@��@��@�  @���@��#@�O�@�(�@��!@�=q@��u@�E�@���@�;d@�~�@��@�x�@��j@�A�@��;@���@��@�@��T@��P@�1'@�33@��^@��j@�1'@��R@���@��D@��;@���@�A�@�1'@���@�E�@��T@�%@��/@�z�@��@�@~ff@}@}�-@}��@{��@x�u@w�w@v�@u�h@t�@s�m@s@r-@qhs@q7L@r�\@r�@p�@pA�@o��@o�@o
=@n��@nE�@m@m�-@m�@k��@i%@h�@g\)@eV@e�@d�@e�T@dI�@b��@dj@a�7@`  @bn�@_�P@^��@_�w@`bN@_
=@^$�@]�-@]�@\�@\�@[dZ@Z��@ZM�@Y��@X��@W�@W�@V��@U�@U?}@T�j@TZ@S�@R�@Qx�@Q7L@Pr�@O�w@O+@N�+@M�T@M?}@L��@K�m@K33@J��@J�@I&�@H�9@Ihs@Ihs@Hr�@FE�@E��@EV@Dz�@C��@C��@C@B^5@A�^@AG�@@�9@@ �@?�w@?;d@>�+@=�@=�h@=O�@<��@;�F@;o@:��@9��@9%@8�9@8A�@7�P@7;d@6�R@6$�@5�-@5�@4��@49X@3�F@3o@2��@2J@1�7@1%@0�u@0b@/l�@/+@.v�@-�T@-`B@,��@,I�@+ƨ@+C�@*~�@)�^@)7L@(��@( �@'��@'\)@&�y@&V@%��@%p�@$�@$j@#�F@#C�@"��@"-@!�^@!7L@ �@ 1'@�w@|�@ff@$�@p�@��@��@j@j@��@�@��@�@�#@��@x�@�`@A�@b@|�@+@��@$�@�h@/@�/@z�@�
@dZ@33@�!@=q@��@x�@&�@�`@�u@bN@  @�w@+@��@5?@��@p�@�@�j@j@1@�F@t�@C�@
��@
^5@
-@	��@	x�@	7L@�`@�u@Q�@b@�w@;d@�@��@5?@@�T@��@O�@�@��@�D@I�@�F@��@t�@o@�H@��@^5@-@�#@x�@7L@%@ �`@ Ĝ@ �@ A�?��;?�\)?��?�v�?���?�O�?��?�1?��?��H?�~�?�=q?���?���?�r�?�1'111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A�v�A�~�A�~�A�z�A�jA�p�A�~�A�hsA��A�p�A�bNA�I�A�I�A�{A~r�Av5?Aex�Ac�-Ab��A^^5AW+ATv�AQ"�AO/AL-AJ$�AI\)AH�AEAD�ABVAA��A>^5A=�;A=�#A=��A<�jA;`BA:~�A9�A8��A8�DA7�wA7"�A77LA5��A5�A61'A6(�A5C�A5A4�A3G�A/��A-��A,�jA+��A+S�A(JA$~�A"VA�^A33A��A�A��A�FAhsA�yAC�@�\)@��T@�z�@��@�{@�?}@��@�j@�x�@��@�`B@�/@�bN@��@�\)@�@�@��@�bN@���@��m@�\)@�@��@���@�Ĝ@�Q�@�@�A�@��H@�h@�r�@�ȴ@��@�(�@�
=@�@�V@��;@�\@�@�dZ@�M�@��T@�%@�ƨ@�"�@�E�@�1'@ա�@��/@�Q�@ӕ�@���@ѡ�@� �@Ώ\@�X@˾w@�K�@�M�@Ȭ@�C�@Ɨ�@��`@�bN@���@���@�I�@��@���@���@���@��@�n�@�5?@�p�@��@�dZ@��y@�$�@�X@�Ĝ@���@�5?@��@�{@��-@��;@�"�@�^5@��@�;d@�@�(�@��#@�n�@��T@��7@�9X@�@��@�=q@�
=@��@��y@�A�@��@��@�5?@�ff@��@��@��@���@�j@���@��@���@�33@��j@���@���@��#@�b@�O�@�p�@�V@��@�  @��@�Ĝ@��@��@��#@�M�@�$�@�`B@�x�@�%@�r�@���@���@��m@���@���@��@��@�  @���@��#@�O�@�(�@��!@�=q@��u@�E�@���@�;d@�~�@��@�x�@��j@�A�@��;@���@��@�@��T@��P@�1'@�33@��^@��j@�1'@��R@���@��D@��;@���@�A�@�1'@���@�E�@��T@�%@��/@�z�@��@�@~ff@}@}�-@}��@{��@x�u@w�w@v�@u�h@t�@s�m@s@r-@qhs@q7L@r�\@r�@p�@pA�@o��@o�@o
=@n��@nE�@m@m�-@m�@k��@i%@h�@g\)@eV@e�@d�@e�T@dI�@b��@dj@a�7@`  @bn�@_�P@^��@_�w@`bN@_
=@^$�@]�-@]�@\�@\�@[dZ@Z��@ZM�@Y��@X��@W�@W�@V��@U�@U?}@T�j@TZ@S�@R�@Qx�@Q7L@Pr�@O�w@O+@N�+@M�T@M?}@L��@K�m@K33@J��@J�@I&�@H�9@Ihs@Ihs@Hr�@FE�@E��@EV@Dz�@C��@C��@C@B^5@A�^@AG�@@�9@@ �@?�w@?;d@>�+@=�@=�h@=O�@<��@;�F@;o@:��@9��@9%@8�9@8A�@7�P@7;d@6�R@6$�@5�-@5�@4��@49X@3�F@3o@2��@2J@1�7@1%@0�u@0b@/l�@/+@.v�@-�T@-`B@,��@,I�@+ƨ@+C�@*~�@)�^@)7L@(��@( �@'��@'\)@&�y@&V@%��@%p�@$�@$j@#�F@#C�@"��@"-@!�^@!7L@ �@ 1'@�w@|�@ff@$�@p�@��@��@j@j@��@�@��@�@�#@��@x�@�`@A�@b@|�@+@��@$�@�h@/@�/@z�@�
@dZ@33@�!@=q@��@x�@&�@�`@�u@bN@  @�w@+@��@5?@��@p�@�@�j@j@1@�F@t�@C�@
��@
^5@
-@	��@	x�@	7L@�`@�u@Q�@b@�w@;d@�@��@5?@@�T@��@O�@�@��@�D@I�@�F@��@t�@o@�H@��@^5@-@�#@x�@7L@%@ �`@ Ĝ@ �@ A�?��;?�\)?��?�v�?���?�O�?��?�1?��?��H?�~�?�=q?���?���?�r�?�1'111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB ��B ��B ��B ��B ��B ��B ��B ��B ��B ��B ��B ��B �hB �=B �hB ��B>wB:^B:^BF�B,B0!B'�B33BQ�Bp�B|�B~�Bm�Bz�B�hB��B��B��B�B��B�}BƨB��B�B�`B�sB��B	7B\B.BG�Bn�B~�B~�B�DB��B��B��B�hB�oB�jBȴB��B�B�yB�B��B.B^5B��B�FB�wB�HB1B:^B\)B� B��B�}B�B�BA�BbNB�1B�BŢB�B��BJB�B5?BbNB�7B�BÖB��B�ZB�yB�B��B	7B�B$�B2-B=qBE�BO�BW
BZB]/B`BBbNBdZBgmBhsBk�Bl�Bk�Bm�Bo�Bo�Bp�Br�Bx�Bx�Bx�By�Bx�B{�Bz�B}�B}�B}�B|�B|�B}�B|�B{�B~�B~�B� B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�+B�=B�hB�oB�uB��B��B��B��B��B�B��B�wB�FBBB�RB��B�^B�?B�!B�B�FB�FB�wBɺB��B��B�
B��B��B��B��B��B��B�B�B�B�
B��B�
B�/B�NB�HB�B�BBB��B�B�mB�mB�mB�B�B�B�B��B��B��BB%B+BDBVBuBuB�B�B!�B$�B'�B<jB?}B=qB<jB9XB8RB>wB?}BB�BF�BJ�BO�BVBZB_;BgmBv�B�B�%B�1B�JB�VB�{B��B�{B��B��B��B��B�B�B�'B�?B�LB�XB�jB�}B�XB��B��B��B��B��B�B�;B�NB�fB�yB�B�B��B��B��B	B	B	1B	JB	VB	{B	�B	�B	�B	$�B	&�B	&�B	&�B	)�B	)�B	-B	1'B	5?B	<jB	<jB	=qB	C�B	A�B	B�B	I�B	H�B	K�B	Q�B	XB	YB	[#B	]/B	_;B	bNB	e`B	gmB	jB	l�B	m�B	o�B	q�B	r�B	v�B	x�B	z�B	z�B	{�B	}�B	~�B	� B	�B	�B	�B	�%B	�+B	�7B	�=B	�JB	�VB	�bB	�oB	�{B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�B	�B	�'B	�3B	�?B	�FB	�FB	�RB	�XB	�^B	�jB	�qB	�}B	��B	��B	ÖB	ĜB	ŢB	ǮB	ȴB	��B	��B	��B	��B	��B	��B	��B	��B	�B	�
B	�B	�#B	�#B	�/B	�5B	�BB	�BB	�HB	�NB	�TB	�ZB	�`B	�mB	�sB	�B	�B	�B	�B	�B	�B	�B	�B	�B	��B	��B	��B	��B	��B	��B	��B	��B	��B
  B
B
B
B
B
B
%B
1B
	7B

=B
DB
JB
PB
VB
VB
\B
hB
uB
{B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
 �B
!�B
"�B
#�B
$�B
$�B
%�B
%�B
&�B
'�B
(�B
)�B
+B
,B
-B
-B
/B
/B
0!B
1'B
2-B
2-B
33B
49B
49B
5?B
6FB
7LB
7LB
9XB
9XB
:^B
:^B
<jB
<jB
<jB
=qB
>wB
>wB
?}B
@�B
A�B
B�B
B�B
C�B
D�B
D�B
D�B
E�B
F�B
G�B
H�B
I�B
I�B
J�B
K�B
L�B
M�B
N�B
N�B
O�B
P�B
P�B
Q�B
R�B
R�B
S�B
T�B
VB
VB
W
B
XB
YB
YB
ZB
[#B
\)111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B �&B �B �B �B �(B �%B �B �&B �B �B �B �B ��B ��B ��B ÂBXBS�BS�B`GBE�BI�BA�BL�Bk�B�HB��B��B�1B��B�B�/B�OB��BĭBçB�*B�WB�B��B�B(B�B"�B)BG�BatB�iB��B��B�B�fB�\B��B�;B�CB�FB�B�B�B]BsB�BHBx-B��B�LB�~B�TB"CBTzBvLB�+B�$BٳB�B7�B[�B|�B��B�{B�BB_B&�B7)BO�B|�B��BŞB�6B�B�B!B(BkB#�B4LB?�BL�BX,B`]Bj�Bq�Bt�Bw�Bz�B}BB�-B�0B�CB�IB�DB�PB�^B�\B�dB�qB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�(B�,B�4B�?B�aB�QB�[B�xB��BľB�?B�
B�VB�YB�B�LB�#B�B��B��B�B�	B�<B�B�B�B��B��B��B�B�B�B�B��B��B��B��B�B��B��B�B�BSB	hB�B�B�B
kB:B5B7BPB	eBpB�B�B�B�B�B �B!�B&B)(B.HB.JB2bB8�B<�B?�BB�BWCBZYBXLBWEBT1BS+BYOBZZB]kBa�Be�Bj�Bp�Bt�BzB�PB��B��B�B�B�1B�@B�gB�kB�cB��B��B��B��B��B�B�B�0B�;B�GB�\B�pB�JB��B��B��B��B��B�B�2B�EB	_B	rB	
�B	�B	�B	�B	�B	B	
B	#1B	'KB	)XB	/}B	3�B	6�B	9�B	?�B	A�B	A�B	A�B	EB	D�B	HB	L/B	PIB	WsB	WsB	X{B	^�B	\�B	]�B	d�B	c�B	f�B	l�B	sB	t%B	v/B	x=B	zIB	}\B	�mB	�|B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�B	�B	�%B	�-B	�;B	�AB	�LB	�UB	�aB	�lB	�wB	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�B	�B	�+B	�5B	�8B	�FB	�RB	�\B	�eB	�dB	�qB	�sB	�|B	׋B	ؐB	ښB	ۣB	ܨB	޷B	߻B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�B	�B	�(B	�,B	�:B	�EB	�EB	�TB	�WB	�eB	�fB	�lB	�sB	�yB	��B
 �B
�B
�B
�B
�B
�B
�B
	�B

�B
�B
�B
�B
�B
�B
�B
�B
�B
B
B
B
#B
)B
,B
7B
;B
EB
 IB
!QB
#\B
$bB
%gB
&mB
'yB
({B
)B
)�B
*�B
,�B
.�B
/�B
0�B
1�B
2�B
4�B
5�B
5�B
6�B
8�B
9�B
:�B
;�B
<�B
=�B
?B
@
B
@
B
AB
AB
BB
CB
D&B
E,B
F3B
G6B
H>B
H>B
JLB
JIB
KSB
LZB
M]B
M\B
NdB
OkB
OkB
PoB
QwB
R�B
R�B
T�B
T�B
U�B
U�B
W�B
W�B
W�B
X�B
Y�B
Y�B
Z�B
[�B
\�B
]�B
]�B
^�B
_�B
_�B
_�B
`�B
a�B
b�B
c�B
d�B
d�B
e�B
f�B
hB
i	B
jB
jB
kB
lB
lB
m$B
n'B
n'B
o0B
p4B
q:B
q:B
rAB
sFB
tNB
tNB
uTB
v]B
wb111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED=PRES + coefficient (see procedure 3.2 in Argo DMQC manual v3.0)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   PRES_ADJUSTED=PRES + coefficient (see procedure 3.2 in Argo DMQC manual v3.3)                                                                                                                                                                                                                                                                                                                                                                                                                                                   PSAL_ADJUSTED is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r.                                                                                                                                                     ADDITIVE COEFFICIENT FOR PRESSURE ADJUSTMENT IS 0.1 dbar                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ADDITIVE COEFFICIENT FOR PRESSURE ADJUSTMENT IS 0.1 dbar                                                                                                                                                                                                                                                                                                                                                                                                                                                                        r=1.000696, +/- 0.000244985                                                                                                                                                                                                                                     PRES_ADJUSTED is calculated following the 3.2 procedure in the Argo Quality Control Manual version 3.0. No significant pressure drift was detected.Pressure evaluation done on 08-Jan-2018 14:20:56                                                             No approved method for delayed-mode qc on TEMP is available                                                                                                                                                                                                     No adjustment is needed on this parameter because no significant sensor drift has been detected.                                                                                                                                                                PRES_ADJUSTED is calculated following the 3.2 procedure in the Argo Quality Control Manual version 3.3. No significant pressure drift was detected.Pressure evaluation done on 15-Mar-2021 09:45:25                                                             No approved method for delayed-mode qc on TEMP is available                                                                                                                                                                                                     Sensor drift/offset detected. Adjusted salinity to OWC(2020) statistical recommendation with CTD_2019V01(WOD2009+),ARGO_2020V01,BOTTLE_2008V1 as reference database. Mapping scales used are 50/4846/53 (lon) 50/4846/53 (lat).                                 201801081427102018010814271020180108142710202103151046142021031510461420210315104614ME  ARDP    1.0                                                                 20150709000000  CR  RCRD            G�O�G�O�G�O�                ME  ARGQ    1.0                                                                 20150709000000  QCF$RCRD            G�O�G�O�G�O�00004000        ME  ARUP    1.0                                                                 20150709000000  UP  RCRD            G�O�G�O�G�O�                ME  JVFM    1.0                                                                 20150709000000  CR  RCRD            G�O�G�O�G�O�                ME  ARCAOW  1.0 CTD_2017V01(WOD2009+),ARGO_2017V02,BOTTLE_2008V1                20151217000000  CV  DOXY            G�O�G�O�G�O�                ME  ARUP    1.0                                                                 20151217000000  UP  RCRD            G�O�G�O�G�O�                ME  ARDU    1.0                                                                 20170824000000  UP  RCRD            G�O�G�O�G�O�                ME  ARSQ    1.1                                                                 20170830000000  QCCVRCRD            G�O�G�O�G�O�                ME  ARGQ                                                                        20170830000000  CF  PSAL            A���A���@�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            A�  A�  @�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            BDffBDff@�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            B�33B�33@�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            B���B���@�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            B�33B�33@�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            CL�CL�@�                  ME  ARGQ                                                                        20170830000000  CF  PSAL            C*L�C*L�@�                  ME  ARSQ    1.1                                                                 20171221000000  QCCVRCRD            G�O�G�O�G�O�                ME  ARDU    1.0                                                                 20171222000000  UP  RCRD            G�O�G�O�G�O�                ME  ARSQ    1.1                                                                 20180108000000  QCCVRCRD            G�O�G�O�G�O�                ME  ARDU    1.0                                                                 20180109000000  UP  RCRD            G�O�G�O�G�O�                ME  ARDU    1.0                                                                 20180529000000  UP  RCRD            G�O�G�O�G�O�                ME  ARSQ    2.0.                                                                20181218000000  QCCVRCRD            G�O�G�O�G�O�                ME  ARDU    1.0                                                                 20181220000000  UP  RCRD            G�O�G�O�G�O�                ME  ARGQ    1.0                                                                 20210315104614  QCP$RCRD            G�O�G�O�G�O�000FFFCE        ME  ARSQOWC 3.0.CTD_2019V01(WOD2009+),ARGO_2020V01,BOTTLE_2008V1                20210315104614  QCCV                G�O�G�O�G�O�                