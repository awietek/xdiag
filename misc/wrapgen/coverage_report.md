# xdiag Julia wrapper — API inventory & coverage

- headers scanned: 58  (parse errors: 0)
- API records (classes/structs): 47
- free functions/operators: 320
- currently add_type<>'d in julia/src: 24

## Type-classification tally (args + returns)

- wrapped: 770
- primitive: 595
- variant: 104
- template_param: 72
- string: 69
- arma: 50
- ostream: 44
- void: 35
- pointer: 31
- stdvector: 20
- std_composite: 14
- unknown: 5

### Unclassified ('unknown') types — need a classifier rule or override

- `std::initializer_list<xdiag::Op>` ×1
- `std::initializer_list<long long>` ×1
- `std::istream` ×1
- `std::function<double ( xdiag::ProductState )>` ×1
- `std::function<std::complex<double> ( xdiag::ProductState )>` ×1

## Internal-type candidates (scope decision)

_Full member lists so you can decide wrap-vs-exclude type-by-type._

### `xdiag::Coeff` — 12 members  (NOT wrapped)

- `Coeff() = default`
- `explicit Coeff(std::string value)`
- `explicit Coeff(const char *value)`
- `explicit Coeff(double value)`
- `explicit Coeff(complex value)`
- `explicit Coeff(Scalar value)`
- `bool operator==(Coeff const &rhs) const`
- `bool operator!=(Coeff const &rhs) const`
- `bool isscalar() const`
- `bool isstring() const`
- `Scalar scalar() const`
- `std::string string() const`

### `xdiag::Matrix` — 27 members  (NOT wrapped)

- `Matrix() = default`
- `Matrix(arma::mat const &mat)`
- `Matrix(arma::cx_mat const &mat)`
- `bool operator==(Matrix const &rhs) const`
- `bool operator!=(Matrix const &rhs) const`
- `bool operator<(Matrix const &rhs) const`
- `Matrix &operator+=(Matrix const &rhs)`
- `Matrix &operator-=(Matrix const &rhs)`
- `Matrix &operator*=(Scalar const &rhs)`
- `Matrix &operator/=(Scalar const &rhs)`
- `Matrix operator-() const`
- `Matrix operator+(Matrix const &b) const`
- `Matrix operator-(Matrix const &b) const`
- `Matrix operator*(Scalar const &rhs) const`
- `Matrix operator/(Scalar const &rhs) const`
- `Matrix operator*(Matrix const &rhs) const`
- `Vector operator*(Vector const &rhs) const`
- `template <typename T> bool is() const`
- `template <typename T> T as() const`
- `int64_t n_rows() const`
- `int64_t n_cols() const`
- `bool isreal() const`
- `arma::mat real() const`
- `arma::mat imag() const`
- `Matrix hc() const`
- `Matrix to_real(double tol = 1e-12) const`
- `bool isapprox(Matrix const &y, double rtol = 1e-12, double atol = 1e-12) const`

### `xdiag::Monomial` — 18 members  (NOT wrapped)

- `Monomial() = default`
- `Monomial(Op const &op)`
- `Monomial(std::initializer_list<Op> ops)`
- `explicit Monomial(std::vector<Op> const &ops)`
- `int64_t size() const noexcept`
- `bool empty() const noexcept`
- `Op const &operator[](int64_t idx) const`
- `std::vector<Op> const &ops() const noexcept`
- `iterator_t begin() const noexcept`
- `iterator_t end() const noexcept`
- `Monomial operator*(Op const &rhs) const`
- `Monomial operator*(Monomial const &rhs) const`
- `Monomial &operator*=(Op const &rhs)`
- `Monomial &operator*=(Monomial const &rhs)`
- `bool operator==(Monomial const &rhs) const noexcept`
- `bool operator!=(Monomial const &rhs) const noexcept`
- `bool operator<(Monomial const &rhs) const noexcept`
- `bool isreal() const`

### `xdiag::Scalar` — 23 members  (NOT wrapped)

- `Scalar() = default`
- `Scalar(double value)`
- `Scalar(complex value)`
- `Scalar(int64_t value)`
- `Scalar(int value)`
- `bool operator==(Scalar const &rhs) const`
- `bool operator!=(Scalar const &rhs) const`
- `Scalar &operator+=(Scalar const &rhs)`
- `Scalar &operator-=(Scalar const &rhs)`
- `Scalar &operator*=(Scalar const &rhs)`
- `Scalar &operator/=(Scalar const &rhs)`
- `Scalar operator-() const`
- `Scalar operator+(Scalar const &b) const`
- `Scalar operator-(Scalar const &b) const`
- `Scalar operator*(Scalar const &b) const`
- `Scalar operator/(Scalar const &b) const`
- `bool isreal() const`
- `double real() const`
- `double imag() const`
- `double abs() const`
- `Scalar conj() const`
- `Scalar to_real( double tol = 1e-12) const`
- `bool isapprox(Scalar const &y, double rtol, double atol) const`

### `xdiag::Term` — 2 members  (NOT wrapped)

- `bool operator==(Term const &rhs) const noexcept`
- `bool operator!=(Term const &rhs) const noexcept`

### `xdiag::Vector` — 25 members  (NOT wrapped)

- `Vector() = default`
- `Vector(arma::vec const &vec)`
- `Vector(arma::cx_vec const &vec)`
- `Vector(std::vector<double> const &vec)`
- `Vector(std::vector<complex> const &vec)`
- `bool operator==(Vector const &rhs) const`
- `bool operator!=(Vector const &rhs) const`
- `Vector &operator+=(Vector const &rhs)`
- `Vector &operator-=(Vector const &rhs)`
- `Vector &operator*=(Scalar const &rhs)`
- `Vector &operator/=(Scalar const &rhs)`
- `Vector operator-() const`
- `Vector operator+(Vector const &b) const`
- `Vector operator-(Vector const &b) const`
- `Vector operator*(Scalar const &rhs) const`
- `Vector operator/(Scalar const &rhs) const`
- `template <typename T> bool is() const`
- `template <typename T> T as() const`
- `int64_t size() const`
- `bool isreal() const`
- `arma::vec real() const`
- `arma::vec imag() const`
- `Vector conj() const`
- `Vector to_real(double tol = 1e-12) const`
- `bool isapprox(Vector const &y, double rtol = 1e-12, double atol = 1e-12) const`

## Distributed / MPI records (excluded from serial Julia wrapper)

- `xdiag::ElectronDistributed`  (xdiag/blocks/distributed/electron_distributed.hpp)
- `xdiag::ElectronDistributedIterator`  (xdiag/blocks/distributed/electron_distributed.hpp)
- `xdiag::SpinhalfDistributed`  (xdiag/blocks/distributed/spinhalf_distributed.hpp)
- `xdiag::SpinhalfDistributedIterator`  (xdiag/blocks/distributed/spinhalf_distributed.hpp)
- `xdiag::tJDistributed`  (xdiag/blocks/distributed/tj_distributed.hpp)
- `xdiag::tJDistributedIterator`  (xdiag/blocks/distributed/tj_distributed.hpp)

## Record coverage

### `xdiag::Boson`  [✗ MISSING] — 15 members

- `Boson() = default`
- `Boson(int64_t sites, int64_t d)`
- `Boson(int64_t sites, int64_t d, int64_t number)`
- `Boson(int64_t sites, int64_t d, Representation const &irrep)`
- `Boson(int64_t sites, int64_t d, int64_t number, Representation const &irrep)`
- `int64_t nsites() const` ✓
- `int64_t d() const`
- `int64_t dim() const` ✓
- `int64_t size() const` ✓
- `bool isreal() const` ✓
- `int64_t index(ProductState const &pstate) const` ✓
- `bool operator==(Boson const &rhs) const` ✓
- `bool operator!=(Boson const &rhs) const` ✓
- `iterator_t begin() const`
- `iterator_t end() const`

### `xdiag::BosonIterator`  [✗ MISSING] — 5 members

- `BosonIterator(Boson const *block, int64_t idx)`
- `BosonIterator &operator++()`
- `ProductState operator*() const` ✓
- `bool operator==(BosonIterator const &rhs) const` ✓
- `bool operator!=(BosonIterator const &rhs) const` ✓

### `xdiag::COOMatrix`  [✓ wrapped] — 0 members


### `xdiag::CSCMatrix`  [✓ wrapped] — 0 members


### `xdiag::CSRMatrix`  [✓ wrapped] — 0 members


### `xdiag::Coeff`  [✗ MISSING] — 12 members

- `Coeff() = default`
- `explicit Coeff(std::string value)`
- `explicit Coeff(const char *value)`
- `explicit Coeff(double value)`
- `explicit Coeff(complex value)`
- `explicit Coeff(Scalar value)`
- `bool operator==(Coeff const &rhs) const` ✓
- `bool operator!=(Coeff const &rhs) const` ✓
- `bool isscalar() const`
- `bool isstring() const`
- `Scalar scalar() const`
- `std::string string() const`

### `xdiag::EigsLanczosResult`  [✗ MISSING] — 6 members

- `arma::vec alphas` ✓
- `arma::vec betas` ✓
- `arma::vec eigenvalues` ✓
- `State eigenvectors` ✓
- `int64_t niterations` ✓
- `std::string criterion` ✓

### `xdiag::EigsLobpcgResult`  [✗ MISSING] — 7 members

- `arma::vec eigenvalues` ✓
- `arma::vec residual_norms`
- `State eigenvectors` ✓
- `int64_t niterations` ✓
- `std::string criterion` ✓
- `arma::mat eigenvalue_history`
- `arma::mat residual_norms_history`

### `xdiag::EigvalsLanczosResult`  [✗ MISSING] — 5 members

- `arma::vec alphas` ✓
- `arma::vec betas` ✓
- `arma::vec eigenvalues` ✓
- `int64_t niterations` ✓
- `std::string criterion` ✓

### `xdiag::Electron`  [✓ wrapped] — 15 members

- `Electron() = default`
- `Electron(int64_t nsites)`
- `Electron(int64_t nsites, int64_t nup, int64_t ndn)`
- `Electron(int64_t nsites, Representation const &irrep)`
- `Electron(int64_t nsites, int64_t nup, int64_t ndn, Representation const &irrep)`
- `int64_t nsites() const` ✓
- `constexpr int64_t d() const`
- `int64_t dim() const` ✓
- `int64_t size() const` ✓
- `bool isreal() const` ✓
- `int64_t index(ProductState const &pstate) const` ✓
- `bool operator==(Electron const &rhs) const` ✓
- `bool operator!=(Electron const &rhs) const` ✓
- `iterator_t begin() const`
- `iterator_t end() const`

### `xdiag::ElectronIterator`  [✗ MISSING] — 5 members

- `ElectronIterator(Electron const *block, int64_t idx)`
- `ElectronIterator &operator++()`
- `ProductState operator*() const` ✓
- `bool operator==(ElectronIterator const &rhs) const` ✓
- `bool operator!=(ElectronIterator const &rhs) const` ✓

### `xdiag::Error`  [✗ MISSING] — 1 members

- `const char *what() const noexcept`

### `xdiag::EvolveLanczosInplaceResult`  [✗ MISSING] — 5 members

- `arma::vec alphas` ✓
- `arma::vec betas` ✓
- `arma::vec eigenvalues` ✓
- `int64_t niterations` ✓
- `std::string criterion` ✓

### `xdiag::EvolveLanczosResult`  [✗ MISSING] — 6 members

- `arma::vec alphas` ✓
- `arma::vec betas` ✓
- `arma::vec eigenvalues` ✓
- `int64_t niterations` ✓
- `std::string criterion` ✓
- `State state` ✓

### `xdiag::Fermion`  [✗ MISSING] — 15 members

- `Fermion() = default`
- `Fermion(int64_t nsites)`
- `Fermion(int64_t nsites, int64_t number)`
- `Fermion(int64_t nsites, Representation const &irrep)`
- `Fermion(int64_t nsites, int64_t number, Representation const &irrep)`
- `int64_t nsites() const` ✓
- `constexpr int64_t d() const`
- `int64_t dim() const` ✓
- `int64_t size() const` ✓
- `bool isreal() const` ✓
- `int64_t index(ProductState const &pstate) const` ✓
- `bool operator==(Fermion const &rhs) const` ✓
- `bool operator!=(Fermion const &rhs) const` ✓
- `iterator_t begin() const`
- `iterator_t end() const`

### `xdiag::FermionIterator`  [✗ MISSING] — 5 members

- `FermionIterator(Fermion const *block, int64_t idx)`
- `FermionIterator &operator++()`
- `ProductState operator*() const` ✓
- `bool operator==(FermionIterator const &rhs) const` ✓
- `bool operator!=(FermionIterator const &rhs) const` ✓

### `xdiag::FileToml`  [✓ wrapped] — 7 members

- `FileToml() = default`
- `explicit FileToml(const char *filename)`
- `explicit FileToml(std::string filename)`
- `explicit FileToml(std::istream &is)`
- `io::FileTomlHandler operator[](std::string key)`
- `bool operator==(FileToml const &other) const` ✓
- `bool operator!=(FileToml const &other) const` ✓

### `xdiag::io::FileTomlHandler`  [✗ MISSING] — 0 members


### `xdiag::GPWF`  [✓ wrapped] — 0 members


### `xdiag::Matrix`  [✗ MISSING] — 27 members

- `Matrix() = default`
- `Matrix(arma::mat const &mat)`
- `Matrix(arma::cx_mat const &mat)`
- `bool operator==(Matrix const &rhs) const` ✓
- `bool operator!=(Matrix const &rhs) const` ✓
- `bool operator<(Matrix const &rhs) const`
- `Matrix &operator+=(Matrix const &rhs)`
- `Matrix &operator-=(Matrix const &rhs)`
- `Matrix &operator*=(Scalar const &rhs)`
- `Matrix &operator/=(Scalar const &rhs)`
- `Matrix operator-() const` ✓
- `Matrix operator+(Matrix const &b) const` ✓
- `Matrix operator-(Matrix const &b) const` ✓
- `Matrix operator*(Scalar const &rhs) const` ✓
- `Matrix operator/(Scalar const &rhs) const` ✓
- `Matrix operator*(Matrix const &rhs) const` ✓
- `Vector operator*(Vector const &rhs) const` ✓
- `template <typename T> bool is() const`
- `template <typename T> T as() const`
- `int64_t n_rows() const` ✓
- `int64_t n_cols() const` ✓
- `bool isreal() const` ✓
- `arma::mat real() const` ✓
- `arma::mat imag() const` ✓
- `Matrix hc() const`
- `Matrix to_real(double tol = 1e-12) const`
- `bool isapprox(Matrix const &y, double rtol = 1e-12, double atol = 1e-12) const` ✓

### `xdiag::Monomial`  [✗ MISSING] — 18 members

- `Monomial() = default`
- `Monomial(Op const &op)`
- `Monomial(std::initializer_list<Op> ops)`
- `explicit Monomial(std::vector<Op> const &ops)`
- `int64_t size() const noexcept` ✓
- `bool empty() const noexcept`
- `Op const &operator[](int64_t idx) const`
- `std::vector<Op> const &ops() const noexcept`
- `iterator_t begin() const noexcept`
- `iterator_t end() const noexcept`
- `Monomial operator*(Op const &rhs) const` ✓
- `Monomial operator*(Monomial const &rhs) const` ✓
- `Monomial &operator*=(Op const &rhs)`
- `Monomial &operator*=(Monomial const &rhs)`
- `bool operator==(Monomial const &rhs) const noexcept` ✓
- `bool operator!=(Monomial const &rhs) const noexcept` ✓
- `bool operator<(Monomial const &rhs) const noexcept`
- `bool isreal() const` ✓

### `xdiag::Op`  [✓ wrapped] — 16 members

- `Op() = default`
- `explicit Op(std::string type)`
- `Op(std::string type, int64_t site)`
- `Op(std::string type, std::vector<int64_t> const &sites)`
- `Op(std::string type, int64_t site, arma::mat const &matrix)`
- `Op(std::string type, std::vector<int64_t> const &sites, arma::mat const &matrix)`
- `Op(std::string type, int64_t site, arma::cx_mat const &matrix)`
- `Op(std::string type, std::vector<int64_t> const &sites, arma::cx_mat const &matrix)`
- `std::string type() const` ✓
- `int64_t size() const` ✓
- `int64_t operator[](int64_t idx) const`
- `std::vector<int64_t> const &sites() const` ✓
- `bool isreal() const` ✓
- `bool operator==(const Op &rhs) const` ✓
- `bool operator!=(const Op &rhs) const` ✓
- `bool operator<(const Op &rhs) const noexcept`

### `xdiag::OpSum`  [✓ wrapped] — 48 members

- `OpSum() = default`
- `explicit OpSum(Op const &op)`
- `explicit OpSum(Monomial const &mono)`
- `OpSum(Scalar const &coeff, Op const &op)`
- `OpSum(Scalar const &coeff, Monomial const &mono)`
- `OpSum(Coeff const &coeff, Op const &op)`
- `OpSum(Coeff const &coeff, Monomial const &mono)`
- `OpSum(std::string const &coeff, Op const &op)`
- `OpSum(double coeff, Op const &op)`
- `OpSum(complex coeff, Op const &op)`
- `OpSum(std::string const &coeff, Monomial const &mono)`
- `OpSum(double coeff, Monomial const &mono)`
- `OpSum(complex coeff, Monomial const &mono)`
- `OpSum &operator+=(OpSum const &ops)`
- `OpSum &operator+=(Op const &op)`
- `OpSum &operator+=(Monomial const &mono)`
- `OpSum operator+(OpSum const &ops) const` ✓
- `OpSum operator+(Op const &op) const` ✓
- `OpSum operator+(Monomial const &mono) const` ✓
- `OpSum &operator-=(OpSum const &ops)`
- `OpSum &operator-=(Op const &op)`
- `OpSum &operator-=(Monomial const &mono)`
- `OpSum operator-(OpSum const &ops) const` ✓
- `OpSum operator-(Op const &op) const` ✓
- `OpSum operator-(Monomial const &mono) const` ✓
- `OpSum operator-() const` ✓
- `OpSum &operator*=(double scalar)`
- `OpSum &operator*=(complex scalar)`
- `OpSum &operator/=(double scalar)`
- `OpSum &operator/=(complex scalar)`
- `OpSum operator*(OpSum const &rhs) const` ✓
- `OpSum &operator*=(OpSum const &rhs)`
- `OpSum operator*(Op const &rhs) const` ✓
- `OpSum &operator*=(Op const &rhs)`
- `OpSum operator*(Monomial const &rhs) const` ✓
- `OpSum &operator*=(Monomial const &rhs)`
- `Scalar &operator[](std::string const &name)`
- `Scalar const &operator[](std::string const &name) const`
- `OpSum plain() const` ✓
- `std::vector<Term> const &terms() const noexcept`
- `std::map<std::string, Scalar> const &params() const noexcept`
- `int64_t size() const noexcept` ✓
- `iterator_t begin() const noexcept`
- `iterator_t end() const noexcept`
- `bool operator==(OpSum const &rhs) const` ✓
- `bool operator!=(OpSum const &rhs) const` ✓
- `bool isreal() const` ✓
- `bool empty() const`

### `xdiag::Permutation`  [✓ wrapped] — 14 members

- `Permutation() = default`
- `explicit Permutation(int64_t size)`
- `Permutation(std::initializer_list<int64_t> list)`
- `explicit Permutation(std::vector<int32_t> const &array)`
- `explicit Permutation(std::vector<int64_t> const &array)`
- `explicit Permutation(arma::Col<int64_t> const &array)`
- `Permutation(int64_t const *ptr, int64_t size)`
- `int64_t size() const` ✓
- `int64_t operator[](int64_t i) const`
- `Permutation inv() const` ✓
- `std::vector<int64_t> const &array() const` ✓
- `Permutation &operator*=(Permutation const &rhs)`
- `bool operator==(Permutation const &rhs) const` ✓
- `bool operator!=(Permutation const &rhs) const` ✓

### `xdiag::PermutationGroup`  [✓ wrapped] — 14 members

- `PermutationGroup() = default`
- `explicit PermutationGroup( std::vector<Permutation> const &permutations)`
- `explicit PermutationGroup( arma::Mat<int64_t> const &matrix)`
- `PermutationGroup(int64_t *ptr, int64_t n_permutations, int64_t nsites)`
- `bool operator==(PermutationGroup const &rhs) const` ✓
- `bool operator!=(PermutationGroup const &rhs) const` ✓
- `int64_t size() const` ✓
- `int64_t nsites() const` ✓
- `Permutation operator[](int64_t sym) const`
- `int64_t const *ptr(int64_t sym) const`
- `int64_t inv(int64_t sym) const` ✓
- `int64_t multiply(int64_t s1, int64_t s2) const` ✓
- `iterator begin() const`
- `iterator end() const`

### `xdiag::Representation::PermutationIrrep`  [✗ MISSING] — 2 members

- `bool operator==(PermutationIrrep const &rhs) const` ✓
- `bool operator!=(PermutationIrrep const &rhs) const` ✓

### `xdiag::ProductState`  [✓ wrapped] — 13 members

- `ProductState() = default`
- `explicit ProductState(int64_t nsites)`
- `explicit ProductState(std::vector<int32_t> const &local_states)`
- `explicit ProductState(std::vector<int64_t> const &local_states)`
- `int64_t operator[](int64_t i) const`
- `int64_t &operator[](int64_t i)`
- `int64_t size() const` ✓
- `int64_t nsites() const` ✓
- `void push_back(int64_t l)`
- `iterator_t begin() const`
- `iterator_t end() const`
- `bool operator==(ProductState const &rhs) const` ✓
- `bool operator!=(ProductState const &rhs) const` ✓

### `xdiag::RandomState`  [✓ wrapped] — 1 members

- `RandomState(int64_t seed = 42, bool normalized = true)`

### `xdiag::Representation`  [✓ wrapped] — 9 members

- `Representation() = default`
- `explicit Representation(PermutationGroup const &group)`
- `Representation(PermutationGroup const &group, Vector const &characters)`
- `Representation(std::string type, int64_t charge)`
- `bool operator==(Representation const &rhs) const` ✓
- `bool operator!=(Representation const &rhs) const` ✓
- `bool is_permutation() const`
- `bool is_charge() const`
- `bool isreal() const` ✓

### `xdiag::Scalar`  [✗ MISSING] — 23 members

- `Scalar() = default`
- `Scalar(double value)`
- `Scalar(complex value)`
- `Scalar(int64_t value)`
- `Scalar(int value)`
- `bool operator==(Scalar const &rhs) const` ✓
- `bool operator!=(Scalar const &rhs) const` ✓
- `Scalar &operator+=(Scalar const &rhs)`
- `Scalar &operator-=(Scalar const &rhs)`
- `Scalar &operator*=(Scalar const &rhs)`
- `Scalar &operator/=(Scalar const &rhs)`
- `Scalar operator-() const` ✓
- `Scalar operator+(Scalar const &b) const` ✓
- `Scalar operator-(Scalar const &b) const` ✓
- `Scalar operator*(Scalar const &b) const` ✓
- `Scalar operator/(Scalar const &b) const` ✓
- `bool isreal() const` ✓
- `double real() const` ✓
- `double imag() const` ✓
- `double abs() const`
- `Scalar conj() const`
- `Scalar to_real( double tol = 1e-12) const`
- `bool isapprox(Scalar const &y, double rtol, double atol) const` ✓

### `xdiag::symmetries::SitePermutation`  [✗ MISSING] — 9 members

- `SitePermutation() = default`
- `explicit SitePermutation(PermutationGroup const &group)`
- `int64_t size() const` ✓
- `int64_t nsites() const` ✓
- `PermutationGroup const &group() const`
- `bool operator==(SitePermutation const &rhs) const` ✓
- `bool operator!=(SitePermutation const &rhs) const` ✓
- `template <typename bit_t> bit_t apply(int64_t sym, bit_t const &bits) const`
- `template <typename bit_t, int nbits> bits::BitArray<bit_t, nbits> apply(int64_t sym, bits::BitArray<bit_t, nbits> const &bits) const`

### `xdiag::Spinhalf`  [✓ wrapped] — 15 members

- `Spinhalf() = default`
- `Spinhalf(int64_t nsites)`
- `Spinhalf(int64_t nsites, int64_t nup)`
- `Spinhalf(int64_t nsites, Representation const &irrep, std::string backend = "auto")`
- `Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep, std::string backend = "auto")`
- `int64_t nsites() const` ✓
- `constexpr int64_t d() const`
- `int64_t dim() const` ✓
- `int64_t size() const` ✓
- `bool isreal() const` ✓
- `int64_t index(ProductState const &pstate) const` ✓
- `bool operator==(Spinhalf const &rhs) const` ✓
- `bool operator!=(Spinhalf const &rhs) const` ✓
- `iterator_t begin() const`
- `iterator_t end() const`

### `xdiag::SpinhalfIterator`  [✗ MISSING] — 5 members

- `SpinhalfIterator(Spinhalf const *block, int64_t idx)`
- `SpinhalfIterator &operator++()`
- `ProductState operator*() const` ✓
- `bool operator==(SpinhalfIterator const &rhs) const` ✓
- `bool operator!=(SpinhalfIterator const &rhs) const` ✓

### `xdiag::State`  [✓ wrapped] — 27 members

- `State() = default`
- `explicit State(Block const &block, bool real = true, int64_t ncols = 1)`
- `State(Block const &block, arma::vec const &vector)`
- `State(Block const &block, arma::cx_vec const &vector)`
- `State(Block const &block, arma::mat const &matrix)`
- `State(Block const &block, arma::cx_mat const &matrix)`
- `bool isvalid() const`
- `int64_t nsites() const` ✓
- `bool isreal() const` ✓
- `State real() const` ✓
- `State imag() const` ✓
- `void make_complex()` ✓
- `int64_t dim() const` ✓
- `int64_t size() const` ✓
- `int64_t nrows() const` ✓
- `int64_t ncols() const` ✓
- `State col(int64_t n, bool copy = true) const` ✓
- `arma::vec vector(int64_t n = 0, bool copy = true) const` ✓
- `arma::mat matrix(bool copy = true) const` ✓
- `arma::cx_vec vectorC(int64_t n = 0, bool copy = true) const` ✓
- `arma::cx_mat matrixC(bool copy = true) const` ✓
- `State(Block const &block, double const *ptr, int64_t ncols, int64_t stride = 1)`
- `State(Block const &block, complex const *ptr, int64_t ncols)`
- `double *memptr()` ✓
- `complex *memptrC()`
- `double *colptr(int64_t col)`
- `complex *colptrC(int64_t col)`

### `xdiag::Term`  [✗ MISSING] — 2 members

- `bool operator==(Term const &rhs) const noexcept` ✓
- `bool operator!=(Term const &rhs) const noexcept` ✓

### `xdiag::TimeEvolveExpokitInplaceResult`  [✗ MISSING] — 2 members

- `double error` ✓
- `double hump` ✓

### `xdiag::TimeEvolveExpokitResult`  [✗ MISSING] — 3 members

- `double error` ✓
- `double hump` ✓
- `State state` ✓

### `xdiag::Vector`  [✗ MISSING] — 25 members

- `Vector() = default`
- `Vector(arma::vec const &vec)`
- `Vector(arma::cx_vec const &vec)`
- `Vector(std::vector<double> const &vec)`
- `Vector(std::vector<complex> const &vec)`
- `bool operator==(Vector const &rhs) const` ✓
- `bool operator!=(Vector const &rhs) const` ✓
- `Vector &operator+=(Vector const &rhs)`
- `Vector &operator-=(Vector const &rhs)`
- `Vector &operator*=(Scalar const &rhs)`
- `Vector &operator/=(Scalar const &rhs)`
- `Vector operator-() const` ✓
- `Vector operator+(Vector const &b) const` ✓
- `Vector operator-(Vector const &b) const` ✓
- `Vector operator*(Scalar const &rhs) const` ✓
- `Vector operator/(Scalar const &rhs) const` ✓
- `template <typename T> bool is() const`
- `template <typename T> T as() const`
- `int64_t size() const` ✓
- `bool isreal() const` ✓
- `arma::vec real() const` ✓
- `arma::vec imag() const` ✓
- `Vector conj() const`
- `Vector to_real(double tol = 1e-12) const`
- `bool isapprox(Vector const &y, double rtol = 1e-12, double atol = 1e-12) const` ✓

### `xdiag::PermutationGroup::iterator`  [✗ MISSING] — 0 members


### `xdiag::tJ`  [✓ wrapped] — 15 members

- `tJ() = default`
- `tJ(int64_t nsites)`
- `tJ(int64_t nsites, int64_t nup, int64_t ndn)`
- `tJ(int64_t nsites, Representation const &irrep)`
- `tJ(int64_t nsites, int64_t nup, int64_t ndn, Representation const &irrep)`
- `int64_t nsites() const` ✓
- `constexpr int64_t d() const`
- `int64_t dim() const` ✓
- `int64_t size() const` ✓
- `bool isreal() const` ✓
- `int64_t index(ProductState const &pstate) const` ✓
- `bool operator==(tJ const &rhs) const` ✓
- `bool operator!=(tJ const &rhs) const` ✓
- `iterator_t begin() const`
- `iterator_t end() const`

### `xdiag::tJIterator`  [✗ MISSING] — 5 members

- `tJIterator(tJ const *block, int64_t idx)`
- `tJIterator &operator++()`
- `ProductState operator*() const` ✓
- `bool operator==(tJIterator const &rhs) const` ✓
- `bool operator!=(tJIterator const &rhs) const` ✓
