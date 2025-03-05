rm -rf qemu/build
mkdir qemu/build
cd qemu/build
../configure --target-list=aarch64-linux-user
make -j$(nproc)
sudo make install
