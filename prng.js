class sha256Random {

HASHLEN = 256 # Number of bits in a hash output
RECIP_HASHLEN = 2**-HASHLEN

    function constructor(seed) {
        if (seed == null) {
            this.seed = Math.random();
        } else {
            this.seed = seed;
        }
        this.counter = 0;
        this.randbits_remaining = 0
        this.hasher = new jsSHA(seed.toString() + ",0", "ASCII");    
    }
    
    function getstate() {
        return(this.seed);
    }
    
    function setseed(seed) { 
        if (seed == null) {
            this.seed = Math.random();
        } else {
            this.seed = seed;
        }
        this.counter = 0;
        this.rand_bits = null;
        this.hasher = new jsSHA(seed.toString() + ",0", "ASCII");
    }
    
    function jumpahead(k) {
        try {
            this.counter += k;
            this.hasher.update("0".repeat(k));
        } catch e {
            console.log(e.toString());
            throw new Error(e);
        }
    }
    
    function nextRandom() {
        r = str2bigInt(this.hasher.getHash("SHA-256", "HEX"), 16, 0)
        this.counter += 1;
        this.hasher.update("0");
        return(r*RECIP_HASHLEN)
    }
    
    function nextRandom() {
        hash_input = (str(self.baseseed) + "," + str(self.counter)).encode('utf-8')
        # Apply SHA-256, interpreting hex output as hexadecimal integer
        # to yield 256-bit integer (a python "long integer")
        hash_output = int(hashlib.sha256(hash_input).hexdigest(),16)
        self.next()
        return(hash_output)
    }
                              
    function getrandbits(k) {
        if (typeof(randbits) == undefined || randbits is null || randbits.len == 0) {                         # initialize the cache
            self.randbits = self.nextRandom()
            self.randbits_remaining = HASHLEN
        while k > self.randbits_remaining{
            self.randbits = (self.nextRandom() << self.randbits_remaining | self.randbits)  # accounts for leading 0s
            self.randbits_remaining = self.randbits_remaining + HASHLEN
        val = (self.randbits & int(2**k-1))            # harvest least significant k bits
        self.randbits_remaining = self.randbits_remaining - k
        self.randbits = self.randbits >> k                 # discard the k harvested bits
        return val
    }
        
    function randbelow_from_randbits(n) {
        k = int(n-1).bit_length()
        r = self.getrandbits(k)   # 0 <= r < 2**k
        while int(r) >= n:
            r = self.getrandbits(k)
        return int(r)
    }        
        
    function randint(a, b, size):
        if size==None:
            return a + self.randbelow_from_randbits(b-a)
        else:
            return np.reshape(np.array([a + self.randbelow_from_randbits(b-a) for i in np.arange(np.prod(size))]), size)

    
    
    function hashMe() {  // hashes the next sequence number
        hasher = new jsSHA($("#seedValue").val() + "," +
                           $("#samNum").val(), "ASCII");
        $("#nObj").val(parseInt($("#nObj").val()));
        if ($("#nObj").val() == 'NaN') {
            $("#nObj").val(1);
        }
        modInt( str2bigInt(hasher.getHash("SHA-256", "HEX"), 16, 0),
                                      parseInt($("#nObj").val())
                                   )
                           );
        }
    }