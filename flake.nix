{
  description = "llvm-pi: compute pi in hand-written LLVM IR";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      supportedSystems = [ "aarch64-darwin" "x86_64-darwin" "aarch64-linux" "x86_64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
    in {
      devShells = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in {
          default = pkgs.mkShell {
            buildInputs = [
              pkgs.gmp
              pkgs.gmp.dev
              pkgs.mpfr
              pkgs.mpfr.dev
              pkgs.llvmPackages_19.llvm
              pkgs.clang_19
            ];
            shellHook = ''
              export GMP_INCLUDE="${pkgs.gmp.dev}/include"
              export GMP_LIB="${pkgs.gmp}/lib"
              export MPFR_INCLUDE="${pkgs.mpfr.dev}/include"
              export MPFR_LIB="${pkgs.mpfr}/lib"
            '';
          };
        });
    };
}
