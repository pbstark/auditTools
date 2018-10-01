class sha256Random {
    constructor(seed) {
        if (seed == null) {
            this.seed = Math.random();
        } else {
            this.seed = seed;
        }
        this.seq = 0;
        this.hasher = new jsSHA(seed.toString() + ",0", "ASCII");    
    }
    
    getseed() {
        return(this.seed);
    }
    
    setseed(seed) { 
        if (seed == null) {
            this.seed = Math.random();
        } else {
            this.seed = seed;
        }
        this.seq = 0;
        this.hasher = new jsSHA(seed.toString() + ",0", "ASCII");
    }
    
    jumpahead(k) {
        try {
            this.seq += k;
            this.hasher.update("0".repeat(k));
        } catch e {
            console.log(e.toString());
        }
    }
    
    nextRandom() {
        r = str2bigInt(this.hasher.getHash("SHA-256", "HEX"), 16, 0)
        this.seq += 1;
        this.hasher.update("0");
        return(r/2**256)
    }
    
    function hashMe() {  // hashes the next sequence number
        hasher = new jsSHA($("#seedValue").val() + "," +
                           $("#samNum").val(), "ASCII");
        $("#nObj").val(parseInt($("#nObj").val()));
        if ($("#nObj").val() == 'NaN') {
            $("#nObj").val(1);
        }
        try {
            sample[0].push($("#samNum").val());
            sample[1].push(hasher.getHash("SHA-256", "HEX"));
            sample[2].push(1 +
                             modInt( str2bigInt(hasher.getHash("SHA-256", "HEX"), 16, 0),
                                      parseInt($("#nObj").val())
                                   )
                           );
            writeList();
            $("#sortedList").val(sample[2].slice().sort(numberLessThan).join(','));
            $("#ballotList").val($("#sortedList").val());
            var deDupeList = sortMultiple(sample[2], numberLessThan);
            $("#sortedDedupeList").val(deDupeList[0].join(','));
            if (vMinMax(deDupeList[1])[1] > 1) {
                 $("#duplicates").val('Ballot, multiplicity\n' + arrayToString(findRepeats(deDupeList)));
            }
        } catch(e) {
        }
    }