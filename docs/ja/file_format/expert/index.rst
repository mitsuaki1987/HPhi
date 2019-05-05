.. highlight:: none

.. _Ch:HowToExpert:

エキスパートモード用入力ファイル
--------------------------------

HPhi++ のエキスパートモードで使用する入力ファイル(\*def)に関して説明します。
入力ファイルの種別は以下の4つで分類されます。

(1) List:
    
    キーワード指定なし: 使用するinput fileの名前のリストを書きます。なお、ファイル名は任意に指定することができます。

(2) Basic parameters:
    
    **CalcMod**: 計算モードを指定するパラメーターを設定します。
    
    **ModPara**: 計算時に必要な基本的なパラメーター(サイトの数、電子数、Lanczosステップを何回やるかなど)を設定します。

    **LocSpin**: 局在スピンの位置を設定します(近藤模型でのみ利用)。

(3) Hamiltonian:
    
    HPhi++ のHamiltonianを電子系の表式により指定します。
    具体的には以下のファイルで指定されます。
    
    | **Trans**:
      :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`\ で表される一体項を指定します。
    | **InterAll**:
      :math:`c_ {i \sigma_1}^{\dagger}c_{j\sigma_2}c_{k \sigma_3}^{\dagger}c_{l \sigma_4}`\ で表される一般二体相互作用を指定します。
	    
    | なお、使用頻度の高い相互作用に関しては下記のキーワードで指定することも可能です。
    | **CoulombIntra**:
      :math:`n_ {i \uparrow}n_{i \downarrow}`\ で表される相互作用を指定します(\ :math:`n_{i \sigma}=c_{i\sigma}^{\dagger}c_{i\sigma}`)。
    | **CoulombInter**:
      :math:`n_ {i}n_{j}`\ で表される相互作用を指定します(\ :math:`n_i=n_{i\uparrow}+n_{i\downarrow}`)。
    | **Hund**:
      :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`\ で表される相互作用を指定します。
    | **PairHop**:
      :math:`c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}`\ で表される相互作用を指定します。
    | **Exchange**:
      :math:`S_i^+ S_j^-`\ で表される相互作用を指定します。  
    | **Ising**: :math:`S_i^z S_j^z`\ で表される相互作用を指定します。  
    | **PairLift**:
      :math:`c_ {i \uparrow}^{\dagger}c_{i\downarrow}c_{j \uparrow}^{\dagger}c_{j \downarrow}`\ で表される相互作用を指定します。

(4) Output:
    
    | **OneBodyG** :出力する一体Green関数を指定します。
      :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`\ が出力されます。
    | **TwoBodyG** :出力する二体Green関数を指定します。
      :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}c^{\dagger}_{k \sigma_3}c_{l\sigma_4}\rangle`

    が出力されます。

.. toctree::
   :maxdepth: 1

   namelist
   calcmod
   modpara
   locspin
   trans
   inter_all
   coulomb_intra
   coulomb_inter
   hund
   pairhop
   exchange
   ising
   pairlift
   green1
   green2
   excitation1
   excitation2
   spectrum_vec
   te1
   te2
