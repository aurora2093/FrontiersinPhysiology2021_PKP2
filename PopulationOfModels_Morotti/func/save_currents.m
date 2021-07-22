function save_currents(i,POM, foldername,PKP2flag,experimflag,BARSflag,CabFac)

    global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store
    global gates Jserca IKs_store Jleak ICFTR Incx Jrel
    global I_ss_store dVm_store Ipca_store I_NaK_store I_Nabk_store I_kr_store
    global I_kur1_store I_kur2_store
    tStep = 1;
    tArray = zeros(1,1e6);
    I_Ca_store=zeros(1,1e6);
    I_to_store=zeros(3,1e6);
    I_Na_store = zeros(1,1e6);
    I_K1_store = zeros(1,1e6);
    ibar_store=zeros(1,1e6);
    gates = zeros(2,1e6);
    Jserca = zeros(1,1e6);
    IKs_store = zeros(1,1e6);
    Jleak = zeros(1e6,2);
    Jrel  = zeros(1,1e6);
    ICFTR = zeros(1,1e6);
    Incx = zeros(1,1e6);
    I_kur1_store = zeros(1,1e6);
    I_kur2_store = zeros(1,1e6);
    I_ss_store = zeros(1,1e6);
    dVm_store = zeros(1,1e6);
    Ipca_store = zeros(1,1e6);
    I_NaK_store = zeros(1,1e6);
    I_Nabk_store = zeros(1,1e6);
    I_kr_store = zeros(1,1e6);
    
    [X0,p] = initCondsMorotti(PKP2flag,BARSflag);
    X0 = POM(i).yfin';
    params = POM(i).params;
    CL = POM(i).CL;
    S = getJacobian;
    options = odeset('RelTol',1e-5,'MaxStep',2, 'JPattern',S);
    
    mod = @myMorotti;
    [time, y]=ode15s(mod,[0 CL],X0,options,p,params,PKP2flag,experimflag,CabFac);
    
    tArray = tArray(1:tStep);
    Ica = I_Ca_store(1:tStep);
    Ito = I_to_store(1,1:tStep);
    Itof = I_to_store(2,1:tStep);
    Itos = I_to_store(3,1:tStep);
    INa = I_Na_store(1:tStep);
    IK1 = I_K1_store(1:tStep);
    s1 = gates(1,1:tStep);
    k1 = gates(2,1:tStep);
    Jserca = Jserca(1:tStep);
    Iks = IKs_store(1:tStep);
    Jleak = Jleak(1:tStep,:);
    Jrel = Jrel(1:tStep);
    ICFTR = ICFTR(1:tStep);
    Incx = Incx(1:tStep);
    Ikur1 = I_kur1_store(1:tStep);
    Ikur2 = I_kur2_store(1:tStep);
    Iss = I_ss_store(1:tStep);
    dVm = dVm_store(1:tStep);
    Ipca = Ipca_store(1:tStep);
    INaK = I_NaK_store(1:tStep);
    INabk = I_Nabk_store(1:tStep);
    Ikr = I_kr_store(1:tStep);
    
    cd(foldername)
    fname = ['PoM' num2str(i) '_currents.mat'];
    save(fname,'time','tArray','Ica','Jserca','Jrel','y','Jleak','Incx','Ito','Itof','Itos','INa','IK1','s1',...
        'k1','Iks','ICFTR','Ikur1','Ikur2','Iss','dVm','Ipca','INaK','INabk','Ikr','y')
    cd F:\Documents\BME\BME_CircAdapt\PKP2_project_cell\SIMULATIONS\POM-master
    
end 