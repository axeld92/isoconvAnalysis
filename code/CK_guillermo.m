function CK_guillermo(u,v)
%% POD-Koopman script algorithm
%
%supposons que dans la matrix u nous avons les donnés (observations) de la component u de la vitesse;
%chaque colonne correspond a chaque temps
%chaque position dans chaque colonne correspond a chaque position dans la
%grille que nous considerons
%le meme pour v 
%generalement,pour la POD, nous supposons que la moyenne de u es 0; le meme
%pour v
%
%M=nombre de observations
%h1=delta x; h2=delta y, du grille
%dt=delta t des observations
%u et v sont NxM, N=nombre de points dans la grille
%%
%c est la matrix de correlation definié par le method des
%snapshots(Sirovich)

%pour la POD

c=u'*u+v'*v;
c=c*h1*h2/M;

[U,S,W]=svd(c,0);

%l=valeurs propres de la POD
l=diag(S);

%nous pouvons choisir combien de modes de POD nous voulons considerer

% traza=sum(l);
% el=cumsum(l)/traza;
% 
% s=find(el>0.9);s=s(1);


at=sqrt(M)*W*sqrt(S);%modes temporelles

%chaque colonne de at correspond a chaque mode temporelle

%components u et v de chaque mode spatial dans chaque colonne

fu=u*at*inv(at'*at);
fv=v*at*inv(at'*at);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pour la DMD

%DMD  si nous considerons les modes temporelles comme observables 



k=at';
[m1,n]=size(k);

k=at';n=M-1;

[Comp residu]=Comp_Koopman(k,n);%avec la matrix companion
[EigVec EigVal]=eig(Comp);
Koop_val=diag(EigVal);

Koop_mod=k(:,1:n)*EigVec;
Koop_init=(Koop_mod'*Koop_mod)\Koop_mod'*k(:,1);
Ki=Koop_init;Kv=Koop_val;Km=Koop_mod;
%vecteur initial, valeurs et modes de Koopman


freq=imag(log(Kv./norm(Kv)))/(2*pi*dt);

for i=1:M;mr(i)=sum(abs(Km(i,:).*Km(i,:)));end
h=l./mr';
E=abs(Km').*abs(Km')*h;%energie

pos=find(freq>0);
loglog(freq(pos),E(pos)/sum(E),'.')


%upk, vpk components u et v des modes POD-Koopman

upk=fu*Km;
vpk=fv*Km;

%%
%%%%%%%%%%%%%%%%%%%%%%%%

%autre maniere

for i=1:n;Va(:,i)=Kv.^(i-1);end;% de Vandermonde
V=k(:,1:n)*inv(Va);%V=modes de Koopman

upk=fu*V;
vpk=fv*V;


for i=1:M;mr(i)=sum(abs(V(i,:).*V(i,:)));end
h=l./mr';
E=abs(V').*abs(V')*h;

pos=find(freq>0);
loglog(freq(pos),E(pos)/sum(E),'.')%le meme resultat
%%%%%%%%%%%%%%%%%%%%




