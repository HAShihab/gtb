(function e$$0(u, s, p) {
    function n(h, r) {
        if (!s[h]) {
            if (!u[h]) {
                var z = "function" == typeof require && require;
                if (!r && z) return z(h, !0);
                if (m) return m(h, !0);
                throw Error("Cannot find module '" + h + "'");
            }
            z = s[h] = {
                exports: {}
            };
            u[h][0].call(z.exports, function(a) {
                var b = u[h][1][a];
                return n(b ? b : a)
            }, z, z.exports, e$$0, u, s, p)
        }
        return s[h].exports
    }
    for (var m = "function" == typeof require && require, h = 0; h < p.length; h++) n(p[h]);
    return n
})({
    1: [function(e, u, s) {
        function p() {}

        function n(d, g, t, q, r) {
            function e(c) {
                if (!c) return q(null, "Couldn't access BAM");
                c = b(c, c.byteLength);
                c = new Uint8Array(c);
                var a = h(c, 0);
                if (a != l) return q(null, "Not a BAM file, magic=0x" + a.toString(16));
                for (var d = h(c, 4), f = "", a = 0; a < d; ++a) f += String.fromCharCode(c[a + 8]);
                f = h(c, d + 8);
                d += 12;
                y.chrToIndex = {};
                y.indexToChr = [];
                for (a = 0; a < f; ++a) {
                    for (var x = h(c, d), g = "", A = 0; A < x - 1; ++A) g += String.fromCharCode(c[d + 4 + A]);
                    h(c, d + x + 4);
                    y.chrToIndex[g] = a;
                    0 == g.indexOf("chr") ? y.chrToIndex[g.substring(3)] = a : y.chrToIndex["chr" + g] = a;
                    y.indexToChr.push(g);
                    d = d + 8 + x
                }
                if (y.indices) return q(y)
            }
            var y = new p;
            y.data = d;
            y.bai =
                g;
            y.indexChunks = t;
            var H = y.indexChunks ? y.indexChunks.minBlockIndex : 1E9;
            if (y.indexChunks) {
                g = y.indexChunks.chunks;
                y.indices = [];
                for (var F = 0; F < g.length; F++) y.indices[F] = null;
                y.data.slice(0, H).fetch(e)
            } else y.bai.fetch(function(c) {
                var b, g, l;
                if (c) {
                    var x = new Uint8Array(c),
                        C = h(x, 0);
                    if (C != f) c = q(null, "Not a BAI file, magic=0x" + C.toString(16));
                    else {
                        C = h(x, 4);
                        y.indices = [];
                        for (var A = 8, J = 0; J < C; ++J) {
                            var k = A;
                            b = x;
                            var E = l = k;
                            g = h(b, E);
                            for (var E = E + 4, P = 0; P < g; ++P) {
                                h(b, E);
                                var K = h(b, E + 4),
                                    E = E + (8 + 16 * K)
                            }
                            for (var P = h(b, E), E = E + 4,
                                    K = 1E9, F = E, m = 0; m < P; ++m) {
                                var z = a(b, F),
                                    F = F + 8;
                                if (z) {
                                    b = z.block;
                                    0 < z.offset && (b += 65536);
                                    b < K && (K = b);
                                    break
                                }
                            }
                            E += 8 * P;
                            b = K;
                            l = E - l;
                            A += l;
                            H = Math.min(b, H);
                            0 < g && (y.indices[J] = new Uint8Array(c, k, A - k))
                        }
                        c = !0
                    }
                } else c = "Couldn't access BAI";
                !0 !== c ? y.bai.url && "undefined" === typeof r ? (y.bai.url = y.data.url.replace(/.bam$/, ".bai"), n(d, y.bai, t, q, !0)) : q(null, c) : y.data.slice(0, H).fetch(e)
            })
        }

        function m() {}
        if ("undefined" !== typeof e) {
            e("./spans");
            s = e("./bin");
            var h = s.readInt,
                v = s.readShort,
                r = s.readByte,
                z = s.readFloat;
            e = e("./lh3utils");
            var a = e.readVob,
                b = e.unbgzf,
                q = e.reg2bins,
                g = e.Chunk
        }
        var l = 21840194,
            f = 21578050;
        e = {
            MULTIPLE_SEGMENTS: 1,
            ALL_SEGMENTS_ALIGN: 2,
            SEGMENT_UNMAPPED: 4,
            NEXT_SEGMENT_UNMAPPED: 8,
            REVERSE_COMPLEMENT: 16,
            NEXT_REVERSE_COMPLEMENT: 32,
            FIRST_SEGMENT: 64,
            LAST_SEGMENT: 128,
            SECONDARY_ALIGNMENT: 256,
            QC_FAIL: 512,
            DUPLICATE: 1024,
            SUPPLEMENTARY: 2048
        };
        p.prototype.blocksForRange = function(d, f, b) {
            var t = this.indices[d];
            if (!t) return [];
            d = q(f, b);
            for (var l = [], r = 0; r < d.length; ++r) l[d[r]] = !0;
            d = [];
            for (var y = [], r = h(t, 0), H = 4, F = 0; F < r; ++F) {
                var c = h(t,
                        H),
                    w = h(t, H + 4),
                    H = H + 8;
                if (l[c])
                    for (var B = 0; B < w; ++B) {
                        var I = a(t, H),
                            x = a(t, H + 8);
                        (4681 > c ? y : d).push(new g(I, x));
                        H += 16
                    } else H += 16 * w
            }
            r = h(t, H);
            l = null;
            f = Math.min(f >> 14, r - 1);
            b = Math.min(b >> 14, r - 1);
            for (r = f; r <= b; ++r)(f = a(t, H + 4 + 8 * r)) && (!l || f.block < l.block || f.offset < l.offset) && (l = f);
            t = [];
            if (null != l)
                for (r = 0; r < y.length; ++r) b = y[r], b.maxv.block >= l.block && b.maxv.offset >= l.offset && t.push(b);
            y = t;
            t = [];
            for (r = 0; r < y.length; ++r) t.push(y[r]);
            for (r = 0; r < d.length; ++r) t.push(d[r]);
            t.sort(function(x, c) {
                var a = x.minv.block - c.minv.block;
                return 0 != a ? a : x.minv.offset - c.minv.offset
            });
            d = [];
            if (0 < t.length) {
                y = t[0];
                for (r = 1; r < t.length; ++r) b = t[r], b.minv.block == y.maxv.block ? y = new g(y.minv, b.maxv) : (d.push(y), y = b);
                d.push(y)
            }
            return d
        };
        p.prototype.fetch = function(a, d, f, g, t) {
            function q() {
                if (B >= F.length) return g(w);
                if (I) {
                    var x = new Uint8Array(I),
                        x = l.readBamRecords(x, F[B].minv.offset, w, d, f, r, t);
                    I = null;
                    ++B;
                    return x ? g(w) : q()
                }
                var c = F[B],
                    x = c.minv.block;
                l.data.slice(x, c.maxv.block + 65536 - x).fetch(function(x) {
                    I = b(x, c.maxv.block - c.minv.block + 1);
                    return q()
                })
            }
            var l = this;
            t = t || {};
            var r = this.chrToIndex[a],
                F;
            if (void 0 === r) F = [];
            else {
                if (null === this.indices[r] && this.indexChunks.chunks[r]) {
                    var c = this.indexChunks.chunks[r];
                    return this.bai.slice(c[0], c[1]).fetch(function(x) {
                        x = new Uint8Array(x);
                        this.indices[r] = x;
                        return this.fetch(a, d, f, g, t)
                    }.bind(this))
                }(F = this.blocksForRange(r, d, f)) || g(null, "Error in index fetch")
            }
            var w = [],
                B = 0,
                I;
            q()
        };
        var d = "=ACxGxxxTxxxxxxN".split(""),
            t = "MIDNSHP=X???????".split("");
        p.prototype.readBamRecords = function(a, f, b, g, q, l, y) {
            for (;;) {
                var H =
                    h(a, f),
                    H = f + H + 4;
                if (H >= a.length) return !1;
                var F = new m,
                    c = h(a, f + 4),
                    w = h(a, f + 8),
                    B = h(a, f + 12),
                    I = (B & 65280) >> 8,
                    x = B & 255,
                    B = h(a, f + 16),
                    C = (B & 4294901760) >> 16,
                    A = B & 65535,
                    B = h(a, f + 20),
                    J = h(a, f + 24),
                    k = h(a, f + 28);
                h(a, f + 32);
                F.segment = this.indexToChr[c];
                F.flag = C;
                F.pos = w;
                F.mq = I;
                y.light && (F.seqLength = B);
                if (!y.light) {
                    0 <= J && (F.nextSegment = this.indexToChr[J], F.nextPos = k);
                    I = "";
                    for (w = 0; w < x - 1; ++w) I += String.fromCharCode(a[f + 36 + w]);
                    F.readName = I;
                    f = f + 36 + x;
                    x = "";
                    for (w = 0; w < A; ++w) I = h(a, f), x = x + (I >> 4) + t[I & 15], f += 4;
                    F.cigar = x;
                    A = "";
                    x = B + 1 >> 1;
                    for (w =
                        0; w < x; ++w) I = a[f + w], A += d[(I & 240) >> 4], A += d[I & 15];
                    f += x;
                    F.seq = A;
                    A = "";
                    for (w = 0; w < B; ++w) A += String.fromCharCode(a[f + w] + 33);
                    f += B;
                    for (F.quals = A; f < H;) {
                        A = String.fromCharCode(a[f], a[f + 1]);
                        x = String.fromCharCode(a[f + 2]);
                        if ("A" == x) x = String.fromCharCode(a[f + 3]), f += 4;
                        else if ("i" == x || "I" == x) x = h(a, f + 3), f += 7;
                        else if ("c" == x || "C" == x) x = a[f + 3], f += 4;
                        else if ("s" == x || "S" == x) x = v(a, f + 3), f += 5;
                        else if ("f" == x) x = z(a, f + 3), f += 7;
                        else if ("Z" == x || "H" == x)
                            for (f += 3, x = ""; w = a[f++], 0 != w;) x += String.fromCharCode(w);
                        else if ("B" == x) {
                            x = String.fromCharCode(a[f +
                                3]);
                            w = h(a, f + 4);
                            if ("i" == x || "I" == x || "f" == x) I = 4, C = "f" == x ? z : h;
                            else if ("s" == x || "S" == x) I = 2, C = v;
                            else if ("c" == x || "C" == x) I = 1, C = r;
                            else throw "Unknown array type " + x;
                            f += 8;
                            x = [];
                            for (J = 0; J < w; ++J) x.push(C(a, f)), f += I
                        } else throw "Unknown type " + x;
                        F[A] = x
                    }
                }
                if (!g || F.pos <= q && F.pos + B >= g) void 0 !== l && c != l || b.push(F);
                if (F.pos > q) return !0;
                f = H
            }
        };
        "undefined" !== typeof u && (u.exports = {
            makeBam: n,
            BAM_MAGIC: l,
            BAI_MAGIC: f,
            BamFlags: e
        })
    }, {
        "./bin": 4,
        "./lh3utils": 24,
        "./spans": 35
    }],
    2: [function(e, u, s) {
        function p(a) {
            this.type = a
        }

        function n(a,
            f) {
            this.parser = a;
            this.sink = f
        }

        function m(a, f) {
            this.parser = a;
            this.sink = f;
            this.wigState = null
        }
        if ("undefined" !== typeof e) {
            var h = e("./spans"),
                v = h.Range,
                r = h.union,
                z = h.intersection,
                h = e("./sourceadapters").registerParserFactory;
            u = e("./das");
            var a = u.DASStylesheet,
                b = u.DASStyle,
                q = u.DASGroup,
                g = e("./utils").shallowCopy
        }
        p.prototype.createSession = function(a) {
            return "wig" == this.type ? new m(this, a) : new n(this, a)
        };
        var l = /([^=]+)=(.+)/,
            f = /\s/,
            d = /^[0-9]+,[0-9]+,[0-9]+/;
        n.prototype.parse = function(a) {
            var b = a.split(f);
            if (!(3 >
                    b.length)) {
                var l = parseInt(b[1]) + 1;
                a = parseInt(b[2]);
                a = {
                    segment: b[0],
                    min: l,
                    max: a
                };
                3 < b.length && "." !== b[3] && (a.label = b[3]);
                4 < b.length && (a.score = parseFloat(b[4]));
                5 < b.length && (a.orientation = b[5]);
                if (8 < b.length) {
                    var h = b[8];
                    d.test(h) && (a.itemRgb = "rgb(" + h + ")")
                }
                if (12 <= b.length) {
                    var h = parseInt(b[6]),
                        e = parseInt(b[7]),
                        m = parseInt(b[9]),
                        n = b[10].split(",").map(function(c) {
                            return parseInt(c)
                        }),
                        y = b[11].split(",").map(function(c) {
                            return parseInt(c)
                        });
                    a.type = "transcript";
                    var H = new q;
                    H.id = b[3];
                    H.type = "transcript";
                    H.notes = [];
                    a.groups = [H];
                    if (12 < b.length) {
                        var F = H = b[12];
                        13 < b.length && (F = b[13]);
                        b = new q;
                        b.id = H;
                        b.label = F;
                        b.type = "gene";
                        a.groups.push(b)
                    }
                    b = null;
                    for (H = 0; H < m; ++H) F = y[H] + l, F = new v(F, F + n[H]), b = b ? r(b, F) : F;
                    y = b.ranges();
                    for (l = 0; l < y.length; ++l) n = y[l], m = g(a), m.min = n.min(), m.max = n.max(), this.sink(m);
                    if (e > h && (h = "+" == a.orientation ? new v(h, e + 3) : new v(h - 3, e), h = z(b, h)))
                        for (a.type = "translation", h = h.ranges(), l = e = 0; l < h.length; ++l) m = l, "-" == a.orientation && (m = h.length - l - 1), n = h[m], m = g(a), m.min = n.min(), m.max = n.max(), a.readframe = e,
                            n = n.max() - n.min(), e = (e + n) % 3, this.sink(m)
                } else this.sink(a)
            }
        };
        n.prototype.flush = function() {};
        m.prototype.parse = function(a) {
            a = a.split(f);
            if ("fixedStep" == a[0]) {
                this.wigState = "fixedStep";
                this.chr = this.pos = this.step = null;
                for (var b = this.span = 1; b < a.length; ++b) {
                    var d = l.exec(a[b]);
                    d && ("chrom" == d[1] ? this.chr = d[2] : "start" == d[1] ? this.pos = parseInt(d[2]) : "step" == d[1] ? this.step = parseInt(d[2]) : "span" == d[1] && (this.span = parseInt(d[2])))
                }
            } else if ("variableStep" == a[0])
                for (this.wigState = "variableStep", this.chr = null, b =
                    this.span = 1; b < a.length; ++b) d = l.exec(a[b]), "chrom" == d[1] ? this.chr = d[2] : "span" == d[1] && (this.span = parseInt(d[2]));
            else this.wigState ? "fixedStep" == this.wigState ? 1 == a.length && (a = parseFloat(a[0]), a = {
                segment: this.chr,
                min: this.pos,
                max: this.pos + this.span - 1,
                score: a
            }, this.pos += this.step, this.sink(a)) : "variableStep" == this.wigState && 2 == a.length && (b = parseInt(a[0]), a = parseFloat(a[1]), a = {
                segment: this.chr,
                min: b,
                max: b + this.span - 1,
                score: a
            }, this.sink(a)) : 4 > a.length || (a = {
                segment: a[0],
                min: parseInt(a[1]) + 1,
                max: parseInt(a[2]),
                score: parseFloat(a[3])
            }, this.sink(a))
        };
        m.prototype.flush = function() {};
        p.prototype.getStyleSheet = function(f) {
            var d = new a;
            if ("wig" == this.type) {
                var g = new b;
                g.glyph = "HISTOGRAM";
                g.BGCOLOR = "blue";
                g.HEIGHT = 30;
                d.pushStyle({
                    type: "default"
                }, null, g)
            } else g = new b, g.glyph = "BOX", g.FGCOLOR = "black", g.BGCOLOR = "blue", g.HEIGHT = 8, g.BUMP = !0, g.LABEL = !0, g.ZINDEX = 20, d.pushStyle({
                type: "default"
            }, null, g), g = new b, g.glyph = "BOX", g.FGCOLOR = "black", g.BGCOLOR = "red", g.HEIGHT = 10, g.BUMP = !0, g.ZINDEX = 20, d.pushStyle({
                    type: "translation"
                },
                null, g), g = new b, g.glyph = "BOX", g.FGCOLOR = "black", g.BGCOLOR = "white", g.HEIGHT = 10, g.ZINDEX = 10, g.BUMP = !0, g.LABEL = !0, d.pushStyle({
                type: "transcript"
            }, null, g), g = new b, g.glyph = "HISTOGRAM", g.COLOR1 = "white", g.COLOR2 = "black", g.HEIGHT = 30, d.pushStyle({
                type: "density"
            }, null, g);
            return f(d)
        };
        h("bed", function(a) {
            return new p(a)
        });
        h("wig", function(a) {
            return new p(a)
        })
    }, {
        "./das": 10,
        "./sourceadapters": 34,
        "./spans": 35,
        "./utils": 48
    }],
    3: [function(e, u, s) {
        function p(a, c) {
            return a[c] + a[c + 1] * L + a[c + 2] * N + a[c + 3] * Q + a[c + 4] * y
        }

        function n() {}

        function m(a, c, f, d) {
            this.bwg = a;
            this.cirTreeOffset = c;
            this.cirTreeLength = f;
            this.isSummary = d
        }

        function h(a, c, f) {
            var d = new n;
            d.data = a;
            d.name = f;
            d.data.slice(0, 512).salted().fetch(function(a) {
                if (!a) return c(null, "Couldn't fetch file");
                var x = new Uint8Array(a),
                    f = new Int16Array(a);
                a = new Int32Array(a);
                var b = x[0] + L * x[1] + N * x[2] + Q * x[3];
                b == t ? d.type = "bigwig" : b == P ? d.type = "bigbed" : b == D || b == G ? c(null, "Currently don't support big-endian BBI files") : c(null, "Not a supported format, magic=0x" + b.toString(16));
                d.version = f[2];
                d.numZoomLevels = f[3];
                d.chromTreeOffset = p(x, 8);
                d.unzoomedDataOffset = p(x, 16);
                d.unzoomedIndexOffset = p(x, 24);
                d.fieldCount = f[16];
                d.definedFieldCount = f[17];
                d.asOffset = p(x, 36);
                d.totalSummaryOffset = p(x, 44);
                d.uncompressBufSize = a[13];
                d.extHeaderOffset = p(x, 56);
                d.zoomLevels = [];
                for (f = 0; f < d.numZoomLevels; ++f) {
                    var b = a[6 * f + 16],
                        g = p(x, 24 * f + 72),
                        k = p(x, 24 * f + 80);
                    d.zoomLevels.push({
                        reduction: b,
                        dataOffset: g,
                        indexOffset: k
                    })
                }
                d.readChromTree(function() {
                    d.getAutoSQL(function(k) {
                        d.schema = k;
                        return c(d)
                    })
                })
            })
        }

        function v(a, c, f,
            d, b) {
            this.bbi = a;
            this.type = c;
            this.fieldCount = f;
            this.offset = d;
            this.field = b
        }
        if ("undefined" !== typeof e) {
            s = e("./spans");
            var r = s.Range,
                z = s.union,
                a = s.intersection;
            s = e("./das");
            var b = s.DASFeature,
                q = s.DASGroup,
                g = e("./utils").shallowCopy,
                l = e("./bin").readInt;
            e = e("jszlib");
            var f = e.inflateBuffer,
                d = e.arrayCopy
        }
        var t = 2291137574,
            D = 654086024,
            P = 2273964779,
            G = 3958540679,
            L = 256,
            N = 65536,
            Q = 16777216,
            y = 4294967296,
            H = /^[0-9]+,[0-9]+,[0-9]+/;
        n.prototype.readChromTree = function(a) {
            var c = this;
            this.chromsToIDs = {};
            this.idsToChroms = {};
            this.maxID = 0;
            var f = this.unzoomedDataOffset,
                f = f + 4 - (f - this.chromTreeOffset & 3);
            this.data.slice(this.chromTreeOffset, f - this.chromTreeOffset).fetch(function(f) {
                var d = new Uint8Array(f),
                    x = new Int16Array(f),
                    b = (new Int32Array(f))[2];
                p(d, 16);
                var A = function(a) {
                    var k = d[a],
                        f = x[a / 2 + 1];
                    a += 4;
                    for (var g = 0; g < f; ++g)
                        if (0 == k) {
                            a += b;
                            var w = p(d, a);
                            a += 8;
                            w -= c.chromTreeOffset;
                            A(w)
                        } else {
                            for (var w = "", l = 0; l < b; ++l) {
                                var q = d[a++];
                                0 != q && (w += String.fromCharCode(q))
                            }
                            l = d[a + 3] << 24 | d[a + 2] << 16 | d[a + 1] << 8 | d[a + 0];
                            a += 8;
                            c.chromsToIDs[w] = l;
                            0 == w.indexOf("chr") && (c.chromsToIDs[w.substr(3)] = l);
                            c.idsToChroms[l] = w;
                            c.maxID = Math.max(c.maxID, l)
                        }
                };
                A(32);
                a(c)
            })
        };
        m.prototype.readWigData = function(a, c, f, d) {
            a = this.bwg.chromsToIDs[a];
            if (void 0 === a) return d([]);
            this.readWigDataById(a, c, f, d)
        };
        m.prototype.readWigDataById = function(a, c, f, d) {
            var b = this;
            if (this.cirHeader) {
                var x = [],
                    g = 0;
                Date.now();
                var A = function(k, x, d, b) {
                        return (0 > a || k == a) && x <= f && d >= c
                    },
                    J = function(c, a) {
                        b.bwg.instrument && console.log("level=" + a + "; offset=" + c + "; time=" + (Date.now() | 0));
                        g += c.length;
                        if (1 == c.length && 48 == c[0] - b.cirTreeOffset && b.cachedCirRoot) E(b.cachedCirRoot, 0, a), --g, 0 == g && b.fetchFeatures(A, x, d);
                        else {
                            for (var f = 4 + 32 * b.cirBlockSize, w, J = 0; J < c.length; ++J) {
                                var l = new r(c[J], c[J] + f);
                                w = w ? z(w, l) : l
                            }
                            f = w.ranges();
                            for (w = 0; w < f.length; ++w) k(c, f[w], a)
                        }
                    },
                    k = function(k, c, a, f) {
                        c.max();
                        c.min();
                        b.bwg.data.slice(c.min(), c.max() - c.min()).fetch(function(f) {
                            for (var w = 0; w < k.length; ++w) c.contains(k[w]) && (E(f, k[w] - c.min(), a), 48 == k[w] - b.cirTreeOffset && 0 == k[w] - c.min() && (b.cachedCirRoot = f), --g, 0 == g && b.fetchFeatures(A,
                                x, d))
                        })
                    },
                    E = function(k, d, b) {
                        var g = new Uint8Array(k),
                            A = new Int16Array(k);
                        k = new Int32Array(k);
                        var E = g[d],
                            A = A[d / 2 + 1];
                        d += 4;
                        if (0 != E)
                            for (E = 0; E < A; ++E) {
                                var l = d / 4,
                                    C = k[l],
                                    q = k[l + 1],
                                    B = k[l + 2],
                                    l = k[l + 3],
                                    t = p(g, d + 16);
                                b = p(g, d + 24);
                                (0 > a || C < a || C == a && q <= f) && (0 > a || B > a || B == a && l >= c) && x.push({
                                    offset: t,
                                    size: b
                                });
                                d += 32
                            } else {
                                for (var I = [], E = 0; E < A; ++E) l = d / 4, C = k[l], q = k[l + 1], B = k[l + 2], l = k[l + 3], t = p(g, d + 16), (0 > a || C < a || C == a && q <= f) && (0 > a || B > a || B == a && l >= c) && I.push(t), d += 24;
                                0 < I.length && J(I, b + 1)
                            }
                    };
                J([b.cirTreeOffset + 48], 1)
            } else this.bwg.data.slice(this.cirTreeOffset,
                48).fetch(function(k) {
                b.cirHeader = k;
                k = new Int32Array(b.cirHeader);
                b.cirBlockSize = k[1];
                b.readWigDataById(a, c, f, d)
            })
        };
        m.prototype.fetchFeatures = function(a, c, g) {
            var l = this;
            c.sort(function(c, x) {
                return (c.offset | 0) - (x.offset | 0)
            });
            if (0 == c.length) g([]);
            else {
                var q = [],
                    x = function(c, x, k, a) {
                        a || (a = {});
                        var f = new b;
                        f._chromId = c;
                        f.segment = l.bwg.idsToChroms[c];
                        f.min = x;
                        f.max = k;
                        f.type = "bigwig";
                        for (var d in a) f[d] = a[d];
                        q.push(f)
                    },
                    C = function() {
                        if (0 == c.length) Date.now(), g(q);
                        else {
                            var b = c[0];
                            if (b.data) l.parseFeatures(b.data,
                                x, a), c.splice(0, 1), C();
                            else {
                                for (var J = b.offset, k = b.size, b = 1; b < c.length && c[b].offset == J + k;) k += c[b].size, ++b;
                                l.bwg.data.slice(J, k).fetch(function(x) {
                                    for (var a = 0, b = 0; a < k;) {
                                        var g = c[b],
                                            A;
                                        0 < l.bwg.uncompressBufSize ? A = f(x, a + 2, g.size - 2) : (A = new Uint8Array(g.size), d(new Uint8Array(x, a, g.size), 0, A, 0, g.size), A = A.buffer);
                                        g.data = A;
                                        a += g.size;
                                        ++b
                                    }
                                    C()
                                })
                            }
                        }
                    };
                C()
            }
        };
        m.prototype.parseFeatures = function(f, c, d) {
            var b = new Uint8Array(f);
            if (this.isSummary)
                for (var l = new Int16Array(f), x = new Int32Array(f), C = new Float32Array(f),
                        b = f.byteLength / 32, l = 0; l < b; ++l) {
                    f = x[8 * l];
                    var A = x[8 * l + 1],
                        J = x[8 * l + 2],
                        k = x[8 * l + 3],
                        E = C[8 * l + 5],
                        t = C[8 * l + 6];
                    d(f, A + 1, J) && (k = {
                        type: "bigwig",
                        score: t / k,
                        maxScore: E
                    }, "bigbed" == this.bwg.type && (k.type = "density"), c(f, A + 1, J, k))
                } else if ("bigwig" == this.bwg.type)
                    if (l = new Int16Array(f), x = new Int32Array(f), C = new Float32Array(f), f = x[0], E = x[1], t = x[3], k = x[4], A = b[20], b = l[11], 3 == A)
                        for (l = 0; l < b; ++l) {
                            var K = C[l + 6],
                                x = E + l * t + 1,
                                A = E + l * t + k;
                            d(f, x, A) && c(f, x, A, {
                                score: K
                            })
                        } else if (2 == A)
                            for (l = 0; l < b; ++l) A = x[2 * l + 6] + 1, J = A + k - 1, K = C[2 * l + 7], d(f, A,
                                J) && c(f, A, J, {
                                score: K
                            });
                        else if (1 == A)
                for (l = 0; l < b; ++l) A = x[3 * l + 6] + 1, J = x[3 * l + 7], K = C[3 * l + 8], A > J && (A = J), d(f, A, J) && c(f, A, J, {
                    score: K
                });
            else console.log("Currently not handling bwgType=" + A);
            else if ("bigbed" == this.bwg.type)
                for (l = 0, k = this.bwg.definedFieldCount, E = this.bwg.schema; l < b.length;) {
                    f = b[l + 3] << 24 | b[l + 2] << 16 | b[l + 1] << 8 | b[l + 0];
                    A = b[l + 7] << 24 | b[l + 6] << 16 | b[l + 5] << 8 | b[l + 4];
                    J = b[l + 11] << 24 | b[l + 10] << 16 | b[l + 9] << 8 | b[l + 8];
                    l += 12;
                    for (K = "";;)
                        if (t = b[l++], 0 != t) K += String.fromCharCode(t);
                        else break;
                    var t = {},
                        h;
                    h = 0 < K.length ?
                        K.split("\t") : [];
                    0 < h.length && 3 < k && (t.label = h[0]);
                    1 < h.length && 4 < k && (K = parseInt(h[1]), isNaN(K) || (t.score = K));
                    2 < h.length && 5 < k && (t.orientation = h[2]);
                    5 < h.length && 8 < k && (K = h[5], H.test(K) && (t.itemRgb = "rgb(" + K + ")"));
                    if (h.length > k - 3 && E)
                        for (K = k - 3; K < h.length; ++K) t[E.fields[K + 3].name] = h[K];
                    if (d(f, A + 1, J, h))
                        if (12 > k) c(f, A + 1, J, t);
                        else {
                            var K = h[3] | 0,
                                D = h[4] | 0,
                                m = h[6] | 0,
                                e = h[7].split(","),
                                y = h[8].split(",");
                            t.exonFrames && (C = t.exonFrames.split(","), t.exonFrames = void 0);
                            t.type = "transcript";
                            var n = new q;
                            for (x in t) n[x] = t[x];
                            n.id = h[0];
                            n.segment = this.bwg.idsToChroms[f];
                            n.min = A + 1;
                            n.max = J;
                            n.notes = [];
                            t.groups = [n];
                            if (9 < h.length) {
                                var p = J = t.geneName || h[9];
                                10 < h.length && (p = h[10]);
                                t.geneName2 && (p = t.geneName2);
                                h = g(n);
                                h.id = J;
                                h.label = p;
                                h.type = "gene";
                                t.groups.push(h)
                            }
                            J = [];
                            for (h = 0; h < m; ++h) n = (y[h] | 0) + A, n = new r(n, n + (e[h] | 0)), J.push(n);
                            e = z(J);
                            m = e.ranges();
                            for (A = 0; A < m.length; ++A) J = m[A], c(f, J.min() + 1, J.max(), t);
                            if (D > K && (A = "+" == t.orientation ? new r(K, D + 3) : new r(K - 3, D), A = a(e, A))) {
                                t.type = "translation";
                                K = A.ranges();
                                for (e = D = 0; K[0].min() > m[e].max();) e++;
                                for (A = 0; A < K.length; ++A) m = A, "-" == t.orientation && (m = K.length - A - 1), J = K[m], t.readframe = D, C && (m = parseInt(C[m + e]), "number" === typeof m && 0 <= m && 2 >= m && (t.readframe = m, t.readframeExplicit = !0)), m = J.max() - J.min(), D = (D + m) % 3, c(f, J.min() + 1, J.max(), t)
                            }
                        }
                } else throw Error("Don't know what to do with " + this.bwg.type);
        };
        m.prototype.getFirstAdjacent = function(a, c, f, d) {
            a = this.bwg.chromsToIDs[a];
            if (void 0 === a) return d([]);
            this.getFirstAdjacentById(a, c, f, d)
        };
        m.prototype.getFirstAdjacentById = function(a, c, f, d) {
            var b = this;
            if (this.cirHeader) {
                var x =
                    null,
                    g = -1,
                    A = -1,
                    l = 0;
                Date.now();
                var k = function(k, c) {
                        l += k.length;
                        for (var x = 4 + 32 * b.cirBlockSize, a, f = 0; f < k.length; ++f) {
                            var d = new r(k[f], k[f] + x);
                            a = a ? z(a, d) : d
                        }
                        x = a.ranges();
                        for (a = 0; a < x.length; ++a) E(k, x[a], c)
                    },
                    E = function(k, g, A, E) {
                        g.max();
                        g.min();
                        b.bwg.data.slice(g.min(), g.max() - g.min()).fetch(function(E) {
                            for (var C = 0; C < k.length; ++C)
                                if (g.contains(k[C]) && (q(E, k[C] - g.min(), A), --l, 0 == l)) {
                                    if (!x) return 0 < f && (0 != a || 0 < c) ? b.getFirstAdjacentById(0, 0, f, d) : 0 > f && (a != b.bwg.maxID || 1E9 > c) ? b.getFirstAdjacentById(b.bwg.maxID,
                                        1E9, f, d) : d([]);
                                    b.fetchFeatures(function(k, x, d, b) {
                                        return 0 > f && (k < a || d < c) || 0 < f && (k > a || x > c)
                                    }, [x], function(k) {
                                        for (var c = null, x = -1, a = -1, b = 0; b < k.length; ++b) {
                                            var g = k[b],
                                                A = g._chromId,
                                                E = g.min,
                                                l = g.max;
                                            if (null == c || 0 > f && (A > x || l > a) || 0 < f && (A < x || E < a)) c = g, a = 0 > f ? l : E, x = A
                                        }
                                        return null != c ? d([c]) : d([])
                                    })
                                }
                        })
                    },
                    q = function(d, E, l) {
                        var J = new Uint8Array(d),
                            q = new Int16Array(d);
                        d = new Int32Array(d);
                        var t = J[E],
                            q = q[E / 2 + 1];
                        E += 4;
                        if (0 != t)
                            for (t = 0; t < q; ++t) {
                                var B = E / 4,
                                    h = d[B],
                                    r = d[B + 1],
                                    m = d[B + 2],
                                    D = d[B + 3],
                                    B = p(J, E + 16);
                                l = p(J, E + 24);
                                (0 > f && (h < a ||
                                    h == a && r <= c) || 0 < f && (m > a || m == a && D >= c)) && !/_random/.exec(b.bwg.idsToChroms[h]) && (null == x || 0 > f && (m > g || m == g && D > A) || 0 < f && (h < g || h == g && r < A)) && (x = {
                                    offset: B,
                                    size: l
                                }, A = 0 > f ? D : r, g = 0 > f ? m : h);
                                E += 32
                            } else {
                                for (var e = J = -1, t = 0; t < q; ++t) {
                                    B = E / 4;
                                    h = d[B];
                                    r = d[B + 1];
                                    m = d[B + 2];
                                    D = d[B + 3];
                                    B = d[B + 4] << 32 | d[B + 5];
                                    if (0 > f && (h < a || h == a && r <= c) && m >= a || 0 < f && (m > a || m == a && D >= c) && h <= a)
                                        if (0 > J || D > e) J = B, e = 0 > f ? D : r;
                                    E += 24
                                }
                                0 <= J && k([J], l + 1)
                            }
                    };
                k([b.cirTreeOffset + 48], 1)
            } else this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(k) {
                b.cirHeader = k;
                k = new Int32Array(b.cirHeader);
                b.cirBlockSize = k[1];
                b.getFirstAdjacentById(a, c, f, d)
            })
        };
        n.prototype.readWigData = function(a, c, f, d) {
            this.getUnzoomedView().readWigData(a, c, f, d)
        };
        n.prototype.getUnzoomedView = function() {
            if (!this.unzoomedView) {
                var a = 4E3;
                this.zoomLevels[0] && (a = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset);
                this.unzoomedView = new m(this, this.unzoomedIndexOffset, a, !1)
            }
            return this.unzoomedView
        };
        n.prototype.getZoomedView = function(a) {
            a = this.zoomLevels[a];
            a.view || (a.view = new m(this, a.indexOffset, 4E3, !0));
            return a.view
        };
        n.prototype._tsFetch = function(a, c, f, d, b) {
            var x = this;
            if (a >= this.zoomLevels.length - 1) {
                if (this.topLevelReductionCache) {
                    for (var g = [], A = this.topLevelReductionCache, l = 0; l < A.length; ++l) A[l]._chromId == c && g.push(A[l]);
                    return b(g)
                }
                this.getZoomedView(this.zoomLevels.length - 1).readWigDataById(-1, 0, 3E8, function(k) {
                    x.topLevelReductionCache = k;
                    return x._tsFetch(a, c, f, d, b)
                })
            } else return (0 > a ? this.getUnzoomedView() : this.getZoomedView(a)).readWigDataById(c, f, d, b)
        };
        n.prototype.thresholdSearch = function(a, c, f, d, b) {
            function x() {
                if (0 ==
                    A.length) return b(null);
                A.sort(function(k, c) {
                    var a = k.zoom - c.zoom;
                    if (0 != a) return a;
                    a = k.chrOrd - c.chrOrd;
                    return 0 != a ? a : k.min - c.min * f
                });
                var k = A.splice(0, 1)[0];
                g._tsFetch(k.zoom, k.chr, k.min, k.max, function(a) {
                    var g = 0 < f ? 0 : 3E8;
                    k.fromRef && (g = c);
                    for (var l = 0; l < a.length; ++l) {
                        var J = a[l],
                            q;
                        q = void 0 != J.maxScore ? J.maxScore : J.score;
                        if (0 < f) {
                            if (q > d)
                                if (0 > k.zoom) {
                                    if (J.min > g) return b(J)
                                } else J.max > g && A.push({
                                    chr: k.chr,
                                    chrOrd: k.chrOrd,
                                    zoom: k.zoom - 2,
                                    min: J.min,
                                    max: J.max,
                                    fromRef: k.fromRef
                                })
                        } else if (q > d)
                            if (0 > k.zoom) {
                                if (J.max <
                                    g) return b(J)
                            } else J.min < g && A.push({
                                chr: k.chr,
                                chrOrd: k.chrOrd,
                                zoom: k.zoom - 2,
                                min: J.min,
                                max: J.max,
                                fromRef: k.fromRef
                            })
                    }
                    x()
                })
            }
            f = 0 > f ? -1 : 1;
            var g = this;
            a = this.chromsToIDs[a];
            for (var A = [{
                    chrOrd: 0,
                    chr: a,
                    zoom: g.zoomLevels.length - 4,
                    min: 0,
                    max: 3E8,
                    fromRef: !0
                }], l = 1; l <= this.maxID + 1; ++l) {
                var k = (a + f * l) % (this.maxID + 1);
                0 > k && (k += this.maxID + 1);
                A.push({
                    chrOrd: l,
                    chr: k,
                    zoom: g.zoomLevels.length - 1,
                    min: 0,
                    max: 3E8
                })
            }
            x()
        };
        n.prototype.getAutoSQL = function(a) {
            if (!this.asOffset) return a(null);
            this.data.slice(this.asOffset, 2048).fetch(function(c) {
                var f =
                    new Uint8Array(c);
                c = "";
                for (var d = 0; d < f.length && 0 != f[d]; ++d) c += String.fromCharCode(f[d]);
                var f = /([\w\[\]]+)\s+(\w+)\s*;\s*("([^"]+)")?\s*/g,
                    b = /(\w+)\s+(\w+)\s+("([^"]+)")?\s+\(\s*/.exec(c);
                if (b) {
                    d = {
                        declType: b[1],
                        name: b[2],
                        comment: b[4],
                        fields: []
                    };
                    c = c.substring(b[0]);
                    for (b = f.exec(c); null != b; b = f.exec(c)) d.fields.push({
                        type: b[1],
                        name: b[2],
                        comment: b[4]
                    });
                    return a(d)
                }
            })
        };
        n.prototype.getExtraIndices = function(a) {
            var c = this;
            if (4 > this.version || 0 == this.extHeaderOffset || "bigbed" != this.type) return a(null);
            this.data.slice(this.extHeaderOffset,
                64).fetch(function(f) {
                if (!f) return a(null, "Couldn't fetch extension header");
                var d = new Uint8Array(f),
                    b = new Int16Array(f);
                new Int32Array(f);
                var x = b[1];
                f = p(d, 4);
                if (0 == x) return a(null);
                c.data.slice(f, 20 * x).fetch(function(f) {
                    if (!f) return a(null, "Couldn't fetch index info");
                    var d = new Uint8Array(f),
                        b = new Int16Array(f);
                    new Int32Array(f);
                    f = [];
                    for (var k = 0; k < x; ++k) {
                        var g = b[10 * k],
                            l = b[10 * k + 1],
                            w = p(d, 20 * k + 4),
                            g = new v(c, g, l, w, b[10 * k + 8]);
                        f.push(g)
                    }
                    a(f)
                })
            })
        };
        v.prototype.lookup = function(a, c) {
            var f = this;
            this.bbi.data.slice(this.offset,
                32).fetch(function(d) {
                function b(k) {
                    f.bbi.data.slice(k, 4 + g * (A + J)).fetch(function(k) {
                        var x = new Uint8Array(k),
                            d = new Uint16Array(k);
                        new Uint32Array(k);
                        d = d[1];
                        k = 4;
                        if (0 == x[0]) {
                            for (var g = null, q = 0; q < d; ++q) {
                                for (var t = "", C = 0; C < A; ++C) {
                                    var h = x[k++];
                                    0 != h && (t += String.fromCharCode(h))
                                }
                                C = p(x, k);
                                k += 8;
                                if (0 > a.localeCompare(t) && g) {
                                    b(g);
                                    return
                                }
                                g = C
                            }
                            b(g)
                        } else {
                            for (q = 0; q < d; ++q) {
                                t = "";
                                for (C = 0; C < A; ++C) h = x[k++], 0 != h && (t += String.fromCharCode(h));
                                if (t == a) return d = p(x, k), x = l(x, k + 8), f.bbi.getUnzoomedView().fetchFeatures(function(k,
                                    c, x, d) {
                                    if (d && d.length > f.field - 3) return d[f.field - 3] == a
                                }, [{
                                    offset: d,
                                    size: x
                                }], c);
                                k += J
                            }
                            return c([])
                        }
                    })
                }
                var x = new Uint8Array(d);
                new Int16Array(d);
                d = new Int32Array(d);
                var g = d[1],
                    A = d[2],
                    J = d[3];
                p(x, 16);
                b(f.offset + 32)
            })
        };
        "undefined" !== typeof u && (u.exports = {
            makeBwg: h,
            BIG_BED_MAGIC: P,
            BIG_WIG_MAGIC: t
        })
    }, {
        "./bin": 4,
        "./das": 10,
        "./spans": 35,
        "./utils": 48,
        jszlib: 64
    }],
    4: [function(e, u, s) {
        function p(a) {
            this.blob = a
        }

        function n(a, d, b, g) {
            g || ("object" === typeof d ? (g = d, d = void 0) : g = {});
            this.url = a;
            this.start = d || 0;
            b && (this.end =
                b);
            this.opts = g
        }

        function m(a) {
            if (!a) return null;
            for (var d = new Uint8Array(a.length), b = 0; b < d.length; ++b) d[b] = a.charCodeAt(b);
            return d.buffer
        }

        function h(a, d) {
            return a[d + 7] << 24 | a[d + 6] << 16 | a[d + 5] << 8 | a[d + 4]
        }

        function v(a, d) {
            return a[d + 3] << 24 | a[d + 2] << 16 | a[d + 1] << 8 | a[d]
        }

        function r(a, d) {
            return a[d + 1] << 8 | a[d]
        }

        function z(a, d) {
            return a[d]
        }

        function a(a, d) {
            return a[d] << 24 | a[d + 1] << 16 | a[d + 2] << 8 | a[d + 3]
        }
        if ("undefined" !== typeof e) var b = e("./utils").shallowCopy,
            q = e("./sha1").b64_sha1;
        p.prototype.slice = function(a, d) {
            var b;
            b = this.blob.slice ? d ? this.blob.slice(a, a + d) : this.blob.slice(a) : d ? this.blob.webkitSlice(a, a + d) : this.blob.webkitSlice(a);
            return new p(b)
        };
        p.prototype.salted = function() {
            return this
        };
        p.prototype.fetch = "undefined" !== typeof FileReader ? function(a) {
            var d = new FileReader;
            d.onloadend = function(b) {
                a(m(d.result))
            };
            d.readAsBinaryString(this.blob)
        } : function(a) {
            var d = new FileReaderSync;
            try {
                var b = d.readAsArrayBuffer(this.blob);
                a(b)
            } catch (g) {
                a(null, g)
            }
        };
        n.prototype.slice = function(a, d) {
            if (0 > a) throw "Bad slice " + a;
            var b =
                this.start,
                g = this.end,
                b = b && a ? b + a : a || b;
            return new n(this.url, b, d && b ? b + d - 1 : g || d - 1, this.opts)
        };
        var g = 0,
            l = 0 <= navigator.userAgent.indexOf("Safari") && 0 > navigator.userAgent.indexOf("Chrome");
        n.prototype.fetchAsText = function(a) {
            var b = new XMLHttpRequest,
                t = this.url;
            (l || this.opts.salt) && 0 > t.indexOf("?") && (t = t + "?salt=" + q("" + Date.now() + "," + ++g));
            b.open("GET", t, !0);
            if (this.end) {
                if (1E8 < this.end - this.start) throw "Monster fetch!";
                b.setRequestHeader("Range", "bytes=" + this.start + "-" + this.end)
            }
            b.onreadystatechange = function() {
                if (4 ==
                    b.readyState) return 200 == b.status || 206 == b.status ? a(b.responseText) : a(null)
            };
            this.opts.credentials && (b.withCredentials = !0);
            b.send("")
        };
        n.prototype.salted = function() {
            var a = b(this.opts);
            a.salt = !0;
            return new n(this.url, this.start, this.end, a)
        };
        n.prototype.fetch = function(a, b, t) {
            var h = this;
            b = b || 1;
            if (3 < b) return a(null);
            var r = new XMLHttpRequest,
                e, n = this.url;
            (l || this.opts.salt) && 0 > n.indexOf("?") && (n = n + "?salt=" + q("" + Date.now() + "," + ++g));
            r.open("GET", n, !0);
            r.overrideMimeType("text/plain; charset=x-user-defined");
            if (this.end) {
                if (1E8 < this.end - this.start) throw "Monster fetch!";
                r.setRequestHeader("Range", "bytes=" + this.start + "-" + this.end);
                e = this.end - this.start + 1
            }
            r.responseType = "arraybuffer";
            r.onreadystatechange = function() {
                if (4 == r.readyState) {
                    if (200 == r.status || 206 == r.status) {
                        if (r.response) {
                            var g = r.response.byteLength;
                            return !e || e == g || t && g == t ? a(r.response) : h.fetch(a, b + 1, g)
                        }
                        if (r.mozResponseArrayBuffer) return a(r.mozResponseArrayBuffer);
                        g = r.responseText;
                        return !e || e == g.length || t && g.length == t ? a(m(r.responseText)) : h.fetch(a,
                            b + 1, g.length)
                    }
                    return h.fetch(a, b + 1)
                }
            };
            this.opts.credentials && (r.withCredentials = !0);
            r.send("")
        };
        (function(a) {
            var b = new ArrayBuffer(8),
                g = new Uint8Array(b),
                l = new Float32Array(b);
            a.readFloat = function(a, b) {
                g[0] = a[b];
                g[1] = a[b + 1];
                g[2] = a[b + 2];
                g[3] = a[b + 3];
                return l[0]
            }
        })(this);
        "undefined" !== typeof u && (u.exports = {
            BlobFetchable: p,
            URLFetchable: n,
            readInt: v,
            readIntBE: a,
            readInt64: h,
            readShort: r,
            readByte: z,
            readFloat: this.readFloat
        })
    }, {
        "./sha1": 33,
        "./utils": 48
    }],
    5: [function(e, u, s) {
        function p(h) {
            var m = "bp";
            1E9 < h ? (h /=
                1E9, m = "Gb") : 1E6 < h ? (h /= 1E6, m = "Mb") : 1E3 < h && (h /= 1E3, m = "kb");
            return "" + Math.round(h) + m
        }
        if ("undefined" !== typeof e) {
            var n = e("./cbrowser").Browser,
                m = e("./utils").makeElement,
                h = e("./numformats").formatLongInt,
                v = e("./zoomslider");
            e("./tier-edit");
            e("./export-config");
            e("./export-ui");
            e("./export-image");
            e("./svg-export");
            e("./session")
        }
        n.prototype.initUI = function(r, e) {
            this.noSourceCSS || ["bootstrap-scoped.css", "dalliance-scoped.css", "font-awesome.min.css"].forEach(function(a) {
                document.head.appendChild(m("link",
                    "", {
                        rel: "stylesheet",
                        href: this.resolveURL("$$css/" + a)
                    }))
            }.bind(this));
            var a = this;
            a.disableDefaultFeaturePopup || this.addFeatureListener(function(c, k, x, b) {
                a.featurePopup(c, k, x, b)
            });
            r.classList.add("dalliance");
            var b = a.toolbar = m("div", null, {
                    className: "btn-toolbar toolbar"
                }),
                q = a.coordSystem.speciesName + " " + a.nameForCoordSystem(a.coordSystem);
            this.setDocumentTitle && (document.title = q + " :: dalliance");
            var g = m("input", "", {
                className: "loc-field"
            });
            a.makeTooltip(g, "Enter a genomic location or gene name");
            var l =
                m("p", "", {
                    className: "loc-status"
                }),
                f = m("a", [m("i", null, {
                    className: "fa fa-search-plus"
                })], {
                    className: "btn"
                }),
                d = new v;
            a.makeTooltip(d, "Highlighted button shows current zoom level, gray button shows inactive zoom level (click or tap SPACE to toggle).");
            var t = m("a", [m("i", null, {
                    className: "fa fa-search-minus"
                })], {
                    className: "btn"
                }),
                D = m("a", [m("i", null, {
                    className: "fa fa-eraser"
                })], {
                    className: "btn"
                }),
                n = m("a", [m("i", null, {
                    className: "fa fa-plus"
                })], {
                    className: "btn"
                }),
                G = m("a", [m("i", null, {
                    className: "fa fa-bookmark"
                })], {
                    className: "btn"
                }),
                L = m("a", [m("i", null, {
                    className: "fa fa-print"
                })], {
                    className: "btn"
                });
            m("a", [m("i", null, {
                className: "fa fa-refresh"
            })], {
                className: "btn"
            });
            var N = m("a", [m("i", null, {
                    className: "fa fa-cogs"
                })], {
                    className: "btn"
                }),
                Q = m("a", [m("i", null, {
                    className: "fa fa-question"
                })], {
                    className: "btn"
                }),
                y = m("a", [m("i", null, {
                    className: "fa fa-road"
                })], {
                    className: "btn"
                });
            a.makeTooltip(y, "Configure currently selected track(s) (E)");
            var H = m("a", [m("i", null, {
                    className: "fa fa-angle-left"
                })], {
                    className: "btn pull-right"
                }, {
                    width: "50px"
                }),
                F = m("a", [m("i", null, {
                    className: "fa fa-angle-right"
                })], {
                    className: "btn pull-right"
                }, {
                    width: "50px"
                }),
                c = m("div", null, {
                    className: "btn-group pull-right"
                });
            this.noTrackAdder || c.appendChild(n);
            this.noTrackEditor || c.appendChild(y);
            this.noExport || c.appendChild(L);
            this.noOptions || c.appendChild(N);
            this.noHelp || c.appendChild(Q);
            this.setUiMode = function(a) {
                this.uiMode = a;
                var k = {
                        help: Q,
                        add: n,
                        opts: N,
                        "export": L,
                        tier: y
                    },
                    c;
                for (c in k) c == a ? k[c].classList.add("active") : k[c].classList.remove("active")
            };
            this.noLeapButtons ||
                b.appendChild(F);
            c.firstChild && b.appendChild(c);
            this.noLeapButtons || b.appendChild(H);
            this.noTitle || b.appendChild(m("div", m("h4", q, {}, {
                margin: "0px"
            }), {
                className: "btn-group title pull-left"
            }));
            this.noLocationField || b.appendChild(m("div", [g, l], {
                className: "btn-group pull-left"
            }));
            this.noClearHighlightsButton || b.appendChild(D);
            this.noZoomSlider || b.appendChild(m("div", [f, m("span", d, {
                className: "btn"
            }), t], {
                className: "btn-group pull-right"
            }));
            this.toolbarBelow ? (r.appendChild(e), r.appendChild(b)) : (r.appendChild(b), r.appendChild(e));
            var w = Math.log(2) / Math.log(10),
                B = Math.log(5) / Math.log(10),
                I = function(c) {
                    var k = (c / a.zoomExpt + Math.log(a.zoomBase)) / Math.log(10);
                    c = k | 0;
                    k -= c;
                    return ((0.01 > k ? c : k <= w + 0.01 ? c + w : k <= B + 0.01 ? c + B : c + 1) * Math.log(10) - Math.log(a.zoomBase)) * a.zoomExpt
                },
                x = function(c) {
                    d.addLabel(c, p(Math.exp(c / a.zoomExpt) * a.zoomBase))
                };
            this.addViewListener(function(c, k, b, A, l) {
                g.value = c + ":" + h(k) + ".." + h(b);
                d.min = l.min | 0;
                d.max = l.max | 0;
                l.isSnapZooming ? (d.value = l.alternate, d.value2 = l.current, d.active = 2) : (d.value = l.current, d.value2 = l.alternate,
                    d.active = 1);
                l.current == l.min ? f.classList.add("disabled") : f.classList.remove("disabled");
                l.current == l.max ? t.classList.add("disabled") : t.classList.remove("disabled");
                d.removeLabels();
                c = l.min;
                l = l.max;
                k = l - c;
                x(I(c));
                x(I(c + 1 * k / 3));
                x(I(c + 2 * k / 3));
                x(I(l));
                a.storeStatus && a.storeViewStatus();
                D.style.display = 0 < a.highlights.length ? "inline-block" : "none"
            });
            this.addTierListener(function() {
                a.storeStatus && a.storeTierStatus()
            });
            g.addEventListener("keydown", function(c) {
                40 == c.keyCode && (c.preventDefault(), c.stopPropagation(),
                    a.setSelectedTier(0));
                if (10 == c.keyCode || 13 == c.keyCode) c.preventDefault(), a.search(g.value, function(a) {
                    l.textContent = a ? "" + a : ""
                })
            }, !1);
            var C;
            n.addEventListener("click", function(c) {
                C && C.displayed ? a.removeAllPopups() : C = a.showTrackAdder(c)
            }, !1);
            a.makeTooltip(n, "Add a new track from the registry or an indexed file. (A)");
            f.addEventListener("click", function(c) {
                c.stopPropagation();
                c.preventDefault();
                a.zoomStep(-10)
            }, !1);
            a.makeTooltip(f, "Zoom in (+)");
            t.addEventListener("click", function(c) {
                c.stopPropagation();
                c.preventDefault();
                a.zoomStep(10)
            }, !1);
            a.makeTooltip(t, "Zoom out (-)");
            d.addEventListener("change", function(c) {
                c = 2 == d.active;
                c != a.isSnapZooming && (a.savedZoom = a.zoomSliderValue - a.zoomMin, a.isSnapZooming = c);
                c = 1 == d.active ? d.value : d.value2;
                a.zoomSliderValue = 1 * c;
                a.zoom(Math.exp(1 * c / a.zoomExpt))
            }, !1);
            G.addEventListener("click", function(a) {
                a.stopPropagation();
                a.preventDefault()
            }, !1);
            a.makeTooltip(G, "Favourite regions");
            L.addEventListener("click", function(c) {
                c.stopPropagation();
                c.preventDefault();
                a.openExportPanel()
            }, !1);
            a.makeTooltip(L, "Export publication-quality SVG. (X)");
            N.addEventListener("click", function(c) {
                c.stopPropagation();
                c.preventDefault();
                a.toggleOptsPopup(c)
            }, !1);
            a.makeTooltip(N, "Configure options.");
            Q.addEventListener("click", function(c) {
                c.stopPropagation();
                c.preventDefault();
                a.toggleHelpPopup(c)
            });
            a.makeTooltip(Q, "Help; Keyboard shortcuts. (H)");
            y.addEventListener("click", function(c) {
                c.stopPropagation();
                c.preventDefault();
                1 == a.selectedTiers.length && a.openTierPanel(a.tiers[a.selectedTiers[0]])
            }, !1);
            H.addEventListener("click", function(c) {
                a.leap(a.reverseKeyScrolling ? -1 : 1, !1)
            }, !1);
            a.makeTooltip(H, function(c) {
                c = a.getSelectedTier();
                var k;
                0 <= c && (k = a.tiers[c]);
                return k && k.featureSource && a.sourceAdapterIsCapable(k.featureSource, "quantLeap") && "number" == typeof k.quantLeapThreshold ? 'Jump to the next region with a score above the threshold in the selected track "' + (k.config.name || k.dasSource.name) + '"" (ctrl+LEFT)' : k && k.featureSource && a.sourceAdapterIsCapable(k.featureSource, "leap") ? 'Jump to the next feature in the selected track "' +
                    (k.config.name || k.dasSource.name) + '" (ctrl+LEFT)' : "Jump left (shift+LEFT)"
            });
            F.addEventListener("click", function(c) {
                a.leap(a.reverseKeyScrolling ? 1 : -1, !1)
            }, !1);
            a.makeTooltip(F, function(c) {
                c = a.getSelectedTier();
                var k;
                0 <= c && (k = a.tiers[c]);
                return k && k.featureSource && a.sourceAdapterIsCapable(k.featureSource, "quantLeap") && "number" == typeof k.quantLeapThreshold ? 'Jump to the next region with a score above the threshold in the selected track "' + (k.config.name || k.dasSource.name) + '"" (ctrl+RIGHT)' : k && k.featureSource &&
                    a.sourceAdapterIsCapable(k.featureSource, "leap") ? 'Jump to the next feature in the selected track "' + (k.config.name || k.dasSource.name) + '" (ctrl+RIGHT)' : "Jump right (shift+RIGHT)"
            });
            a.addTierSelectionListener(function() {
                var c = a.getSelectedTier(),
                    k;
                0 <= c && (k = a.tiers[c]);
                c = !1;
                k && k.featureSource && (a.sourceAdapterIsCapable(k.featureSource, "quantLeap") && "number" == typeof k.quantLeapThreshold ? c = !0 : a.sourceAdapterIsCapable(k.featureSource, "leap") && (c = !0));
                H.firstChild.className = c ? "fa fa-angle-double-left" : "fa fa-angle-left";
                F.firstChild.className = c ? "fa fa-angle-double-right" : "fa fa-angle-right"
            });
            D.addEventListener("click", function(c) {
                a.clearHighlights()
            }, !1);
            a.makeTooltip(D, "Clear highlights (C)");
            a.addTierSelectionWrapListener(function(c) {
                0 > c && (a.setSelectedTier(null), g.focus())
            });
            a.addTierSelectionListener(function(c) {
                "tier" === a.uiMode && (0 == c.length ? (a.hideToolPanel(), a.manipulatingTier = null, a.uiMode = "none") : (c = a.tiers[c[0]], c != a.manipulatingTier && a.openTierPanel(c)))
            });
            var A = function(c) {
                if (65 == c.keyCode || 97 == c.keyCode) c.preventDefault(),
                    c.stopPropagation(), a.showTrackAdder();
                else if (72 == c.keyCode || 104 == c.keyCode) c.stopPropagation(), c.preventDefault(), a.toggleHelpPopup(c);
                else if (69 == c.keyCode || 101 == c.keyCode) c.stopPropagation(), c.preventDefault(), 1 == a.selectedTiers.length && a.openTierPanel(a.tiers[a.selectedTiers[0]]);
                else if (88 == c.keyCode || 120 == c.keyCode) c.stopPropagation(), c.preventDefault(), a.openExportPanel();
                else if (67 == c.keyCode || 99 == c.keyCode) c.stopPropagation(), c.preventDefault(), a.clearHighlights()
            };
            r.addEventListener("focus",
                function(c) {
                    r.addEventListener("keydown", A, !1)
                }, !1);
            r.addEventListener("blur", function(c) {
                r.removeEventListener("keydown", A, !1)
            }, !1);
            r.addEventListener("keydown", function(c) {
                27 === c.keyCode && "none" !== a.uiMode && (c.preventDefault(), c.stopPropagation(), a.setUiMode("none"), a.hideToolPanel(), a.selectedTiers && 0 < a.selectedTiers.length && a.browserHolder.focus())
            }, !1)
        };
        n.prototype.showToolPanel = function(h, e) {
            var a = this;
            this.activeToolPanel && this.activeToolPanel.parentElement.removeChild(this.activeToolPanel);
            var b;
            b = e ? h : m("div", h, {}, {
                overflowY: "auto",
                width: "100%"
            });
            var q = m("div", m("i", null, {
                className: "fa fa-caret-right"
            }), {
                className: "tool-divider"
            });
            q.addEventListener("click", function(b) {
                a.hideToolPanel();
                a.setUiMode("none")
            }, !1);
            this.makeTooltip(q, "Close tool panel (ESC)");
            this.activeToolPanel = m("div", [q, b], {
                className: "tool-holder"
            });
            this.svgHolder.appendChild(this.activeToolPanel);
            this.resizeViewer();
            a = this
        };
        n.prototype.hideToolPanel = function() {
            this.activeToolPanel && this.activeToolPanel.parentElement.removeChild(this.activeToolPanel);
            this.svgHolder.style.width = "100%";
            this.activeToolPanel = null;
            this.resizeViewer()
        };
        n.prototype.toggleHelpPopup = function(h) {
            "help" === this.uiMode ? (this.hideToolPanel(), this.setUiMode("none")) : (h = m("iframe", null, {
                scrolling: "yes",
                seamless: "seamless",
                src: this.resolveURL("$$help/index.html"),
                className: "help-panel"
            }), this.showToolPanel(h, !1), this.setUiMode("help"))
        };
        n.prototype.toggleOptsPopup = function(h) {
            var e = this;
            if ("opts" === this.uiMode) this.hideToolPanel(), this.setUiMode("none");
            else {
                h = m("div", null, {
                    className: "form-horizontal"
                }, {
                    boxSizing: "border-box",
                    MozBoxSizing: "border-box",
                    display: "inline-block",
                    verticalAlign: "top"
                });
                var a = m("table");
                a.cellPadding = 5;
                var b = m("input", "", {
                    type: "checkbox",
                    checked: e.reverseScrolling
                });
                b.addEventListener("change", function(a) {
                    e.reverseScrolling = b.checked;
                    e.storeStatus()
                }, !1);
                a.appendChild(m("tr", [m("td", "Reverse trackpad scrolling", {
                    align: "right"
                }), m("td", b)]));
                var q = m("input", "", {
                    type: "checkbox",
                    checked: e.reverseKeyScrolling
                });
                q.addEventListener("change", function(a) {
                    e.reverseKeyScrolling =
                        q.checked;
                    e.storeStatus()
                }, !1);
                a.appendChild(m("tr", [m("td", "Reverse scrolling buttons and keys", {
                    align: "right"
                }), m("td", q)]));
                var g = m("select");
                g.appendChild(m("option", "Left", {
                    value: "left"
                }));
                g.appendChild(m("option", "Center", {
                    value: "center"
                }));
                g.appendChild(m("option", "Right", {
                    value: "right"
                }));
                g.appendChild(m("option", "None", {
                    value: "none"
                }));
                g.value = e.rulerLocation;
                g.addEventListener("change", function(a) {
                    e.rulerLocation = g.value;
                    e.positionRuler();
                    for (a = 0; a < e.tiers.length; ++a) e.tiers[a].paintQuant();
                    e.storeStatus()
                }, !1);
                a.appendChild(m("tr", [m("td", "Vertical guideline", {
                    align: "right"
                }), m("td", g)]));
                var l = m("input", "", {
                    type: "checkbox",
                    checked: e.singleBaseHighlight
                });
                l.addEventListener("change", function(a) {
                    e.singleBaseHighlight = l.checked;
                    e.positionRuler();
                    e.storeStatus()
                }, !1);
                l.setAttribute("id", "singleBaseHightlightButton");
                a.appendChild(m("tr", [m("td", "Display and highlight current genome location", {
                    align: "right"
                }), m("td", l)]));
                h.appendChild(a);
                a = m("button", "Reset browser", {
                    className: "btn"
                }, {
                    marginLeft: "auto",
                    marginRight: "auto",
                    display: "block"
                });
                a.addEventListener("click", function(a) {
                    e.reset()
                }, !1);
                h.appendChild(a);
                this.showToolPanel(h);
                this.setUiMode("opts")
            }
        }
    }, {
        "./cbrowser": 6,
        "./export-config": 14,
        "./export-image": 15,
        "./export-ui": 16,
        "./numformats": 26,
        "./session": 32,
        "./svg-export": 37,
        "./tier-edit": 43,
        "./utils": 48,
        "./zoomslider": 51
    }],
    6: [function(e, u, s) {
        function p(c, a, b) {
            this.min = a;
            this.max = b;
            this.chr = c
        }

        function n(c) {
            c || (c = {});
            this.prefix = "//www.biodalliance.org/release-0.13/";
            this.sources = [];
            this.tiers = [];
            this.tierGroups = {};
            this.featureListeners = [];
            this.featureHoverListeners = [];
            this.viewListeners = [];
            this.regionSelectListeners = [];
            this.tierListeners = [];
            this.tierSelectionListeners = [];
            this.tierSelectionWrapListeners = [];
            this.cookieKey = "browser";
            this.registry = "//www.dasregistry.org/das/sources";
            this.chains = {};
            this.pageName = "svgHolder";
            this.maxExtra = 2.5;
            this.minExtra = 0.5;
            this.zoomFactor = 1;
            this.maxPixelsPerBase = 10;
            this.origin = 0;
            this.targetQuantRes = 1;
            this.featurePanelWidth = 750;
            this.zoomBase = 100;
            this.zoomExpt = 30;
            this.zoomSliderValue = 100;
            this.entryPoints = null;
            this.currentSeqMax = -1;
            this.highlights = [];
            this.selectedTiers = [1];
            //this.maxViewWidth = 5E5;
            this.maxViewWidth = 1E5;
            this.defaultSubtierMax = 100;
            //this.reverseScrolling = !1;
            this.reverseScrolling = 1;
            this.rulerLocation = "center";
            this.defaultHighlightFill = "red";
            this.defaultHighlightAlpha = 0.3;
            this.singleBaseHighlight = this.exportRuler = this.exportHighlights = !0;
            this.tierBackgroundColors = ["rgb(245,245,245)", "white"];
            this.minTierHeight = 20;
            this.noDefaultLabels = !1;
            this.availableSources = new l;
            this.defaultSources = [];
            this.mappableSources = {};
            this.hubs = [];
            this.hubObjects = [];
            this.sourceCache = new g;
            this.useFetchWorkers = this.retina = !0;
            this.maxWorkers = 2;
            this.workerPath = "$$worker-all.js";
            this.assemblyNameUcsc = this.assemblyNamePrimary = !0;
            this.initListeners = [];
            this.baseColors = c.baseColors ? c.baseColors : {
                A: "green",
                C: "blue",
                G: "orange",
                T: "red",
                "-": "hotpink",
                I: "red"
            };
            if (void 0 !== c.viewStart && "number" !== typeof c.viewStart) throw Error("viewStart must be an integer");
            if (void 0 !== c.viewEnd && "number" !== typeof c.viewEnd) throw Error("viewEnd must be an integer");
            for (var a in c) this[a] = c[a];
            "string" === typeof c.uiPrefix && "string" !== typeof c.prefix && (this.prefix = c.uiPrefix);
            0 > this.prefix.indexOf("//") && 0 === this.prefix.indexOf("/") && (this.prefix = "//" + window.location.hostname + this.prefix);
            0 === this.prefix.indexOf("//") && (c = window.location.protocol, "http:" != c && "https:" != c && (console.log(window.location.protocol), console.log("WARNING: prefix is set to a protocol-relative URL (" + this.prefix + " when loading from a non-HTTP source"), this.prefix = "http:" + this.prefix));
            if (!this.coordSystem) throw Error("Coordinate system must be configured");
            if (void 0 === this.chr || void 0 === this.viewStart || void 0 === this.viewEnd) throw Error("Viewed region (chr:start..end) must be defined");
            var b = this;
            "complete" === document.readyState ? b.realInit() : window.addEventListener("load", function(c) {
                b.realInit()
            }, !1)
        }

        function m(c, a, b) {
            if (b)
                for (var d = 0; d < a.length; ++d) a[d].mapping = b;
            c.set(a)
        }

        function h(c) {
            return c.uri ? c.uri : c.blob ? "file:" + c.blob.name : c.bwgBlob ? "file:" + c.bwgBlob.name : c.bamBlob ? "file:" +
                c.bamBlob.name : c.twoBitBlob ? "file:" + c.twoBitBlob.name : c.bwgURI || c.bamURI || c.jbURI || c.twoBitURI || "http://www.biodalliance.org/magic/no_uri"
        }

        function v(c) {
            return c.stylesheet_uri ? c.stylesheet_uri : "sequence" == c.tier_type || c.twoBitURI || c.twoBitBlob ? "http://www.biodalliance.org/magic/sequence" : h(c)
        }

        function r(c, a) {
            if (h(c) != h(a) || c.mapping != a.mapping || c.tier_type != a.tier_type) return !1;
            if (c.overlay) {
                if (!a.overlay || a.overlay.length != c.overlay.length) return !1;
                for (var b = 0; b < c.overlay.length; ++b)
                    if (!r(c.overlay[b],
                            a.overlay[b])) return !1
            } else if (a.overlay) return !1;
            return !0
        }

        function z(c, a) {
            if (h(c) != h(a) || v(c) != v(a) || c.mapping != a.mapping || c.tier_type != a.tier_type) return !1;
            if (c.overlay) {
                if (!a.overlay || a.overlay.length != c.overlay.length) return !1;
                for (var b = 0; b < c.overlay.length; ++b)
                    if (!z(c.overlay[b], a.overlay[b])) return !1
            } else if (a.overlay) return !1;
            return !0
        }

        function a(c, b, d, f) {
            f = f || [];
            for (var k = c.length - 1; 0 <= k; --k) {
                var g = c[k];
                if (!g.notSelectable && g.min() <= b && g.max() >= b && (!g.minY || !(d < g.minY() || d > g.maxY()))) {
                    g.feature ?
                        f.push(g.feature) : g.group && f.push(g.group);
                    if (g.glyphs) return a(g.glyphs, b, d, f);
                    if (g.glyph) return a([g.glyph], b, d, f);
                    break
                }
            }
            return f
        }

        function b(c, a) {
            var b = this;
            this.tagSeed = 0;
            this.callbacks = {};
            this.browser = c;
            this.worker = a;
            this.worker.onmessage = function(c) {
                var a = b.callbacks[c.data.tag];
                a && (a(c.data.result, c.data.error), delete b.callbacks[c.data.tag])
            }
        }

        function q(c) {
            var a = c.resolveURL(c.workerPath);
            0 == a.indexOf("//") && (a = "https:" == window.location.protocol ? "https:" + a : "http:" + a);
            var d = new Blob(['importScripts("' +
                a + '");'
            ], {
                type: "application/javascript"
            });
            return new F(function(a, k) {
                var f = new Worker(URL.createObjectURL(d));
                f.onmessage = function(k) {
                    "init" === k.data.tag && (console.log("Worker initialized"), a(new b(c, f)))
                };
                f.onerror = function(c) {
                    k(c.message)
                }
            })
        }

        function g() {
            this.sourcesByURI = {}
        }
        if ("undefined" !== typeof e) {
            s = e("./utils");
            var l = s.Observed,
                f = s.makeElement,
                d = s.removeChildren,
                t = s.miniJSONify,
                D = s.shallowCopy,
                P = e("./tier").DasTier,
                G = e("./sha1").hex_sha1,
                L = e("./thub").connectTrackHub,
                N = e("./version");
            s = e("./numformats");
            var Q = s.formatQuantLabel,
                y = s.formatLongInt,
                H = e("./chainset").Chainset,
                F = e("es6-promise").Promise
        }
        n.prototype.resolveURL = function(c) {
            return c.replace("$$", this.prefix)
        };
        n.prototype.realInit = function() {
            var c = this;
            if (this.wasInitialized) console.log("Attemping to call realInit on an already-initialized Dalliance instance");
            else {
                this.wasInitialized = !0;
                var a = navigator.userAgent || "dummy";
                0 <= a.indexOf("Trident") && 0 <= a.indexOf("rv:11") && (this.disablePinning = !0);
                this.defaultChr = this.chr;
                this.defaultStart = this.viewStart;
                this.defaultEnd = this.viewEnd;
                this.defaultSources = [];
                for (a = 0; a < this.sources.length; ++a) {
                    var b = this.sources[a];
                    b && this.defaultSources.push(b)
                }
                this.restoreStatus && (this.statusRestored = this.restoreStatus());
                var g = this;
                this.browserHolderHolder = document.getElementById(this.pageName);
                this.browserHolderHolder.classList.add("dalliance-injection-point");
                this.browserHolder = f("div", null, {
                    className: "dalliance dalliance-root",
                    tabIndex: -1
                });
                this.maxHeight ? this.browserHolder.style.maxHeight = this.maxHeight + "px" : void 0 !=
                    this.maxHeight && (this.browserHolder.style.maxHeight = null);
                d(this.browserHolderHolder);
                this.browserHolderHolder.appendChild(this.browserHolder);
                this.svgHolder = f("div", null, {
                    className: "main-holder"
                });
                this.initUI(this.browserHolder, this.svgHolder);
                this.pinnedTierHolder = f("div", null, {
                    className: "tier-holder tier-holder-pinned"
                });
                this.tierHolder = f("div", this.makeLoader(24), {
                    className: "tier-holder tier-holder-rest"
                });
                this.locSingleBase = f("span", "", {
                    className: "loc-single-base"
                });
                a = f("div", this.locSingleBase, {
                    className: "loc-single-base-holder"
                });
                this.addViewListener(function(a, k, b, d, f, g, l) {
                    k = Math.round((l + g) / 2);
                    c.locSingleBase.appendChild(document.createTextNode(a + ":" + y(k)));
                    c.locSingleBase.removeChild(c.locSingleBase.firstChild)
                });
                this.disablePinning ? this.tierHolderHolder = this.tierHolder : (this.tierHolderHolder = f("div", [a, this.pinnedTierHolder, this.tierHolder], {
                    className: "tier-holder-holder"
                }), this.svgHolder.appendChild(this.tierHolderHolder));
                this.svgHolder.appendChild(this.tierHolderHolder);
                this.bhtmlRoot =
                    f("div");
                this.disablePoweredBy || this.bhtmlRoot.appendChild(f("span", ["Powered by ", f("a", "Biodalliance", {
                    href: "http://www.biodalliance.org/"
                }), " " + N], {
                    className: "powered-by"
                }));
                this.browserHolder.appendChild(this.bhtmlRoot);
                window.addEventListener("resize", function(c) {
                    g.resizeViewer()
                }, !1);
                this.ruler = f("div", null, {
                    className: "guideline"
                });
                this.ruler2 = f("div", null, {
                    className: "single-base-guideline"
                });
                this.tierHolderHolder.appendChild(this.ruler);
                this.tierHolderHolder.appendChild(this.ruler2);
                this.chainConfigs =
                    this.chains || {};
                this.chains = {};
                for (var k in this.chainConfigs) a = this.chainConfigs[k], a instanceof H && console.log('WARNING: Should no longer use "new Chainset" in Biodalliance configurations.'), this.chains[k] = new H(a);
                if (0 < this.maxWorkers) {
                    k = [];
                    for (a = 0; a < this.maxWorkers; ++a) k.push(q(this));
                    k = F.all(k)
                } else k = F.resolve([]);
                this.fetchWorkers = null;
                this.nextWorker = 0;
                k.then(function(c) {
                    console.log("Booted " + c.length + " workers");
                    g.fetchWorkers = c
                }, function(c) {
                    console.log("Failed to boot workers", c)
                }).then(function() {
                    if ("none" !=
                        window.getComputedStyle(g.browserHolderHolder).display && 0 < g.tierHolder.getBoundingClientRect().width) setTimeout(function() {
                        g.realInit2()
                    }, 1);
                    else var c = setInterval(function() {
                        "none" != window.getComputedStyle(g.browserHolderHolder).display && 0 < g.tierHolder.getBoundingClientRect().width && (clearInterval(c), g.realInit2())
                    }, 300)
                })
            }
        };
        n.prototype.realInit2 = function() {
            var c = this;
            d(this.tierHolder);
            d(this.pinnedTierHolder);
            this.featurePanelWidth = this.tierHolder.getBoundingClientRect().width | 0;
            this.scale = this.featurePanelWidth /
                (this.viewEnd - this.viewStart);
            this.zoomMax || (this.zoomMax = this.zoomExpt * Math.log(this.maxViewWidth / this.zoomBase), this.zoomMin = this.zoomExpt * Math.log(this.featurePanelWidth / this.maxPixelsPerBase / this.zoomBase));
            this.zoomSliderValue = this.zoomExpt * Math.log((this.viewEnd - this.viewStart + 1) / this.zoomBase);
            this.tierHolderHolder.addEventListener("mousewheel", function(a) {
                a.stopPropagation();
                a.preventDefault();
                if (a.wheelDeltaX) {
                    var k = a.wheelDeltaX / 5;
                    c.reverseScrolling || (k = -k);
                    c.move(k)
                }
                a.wheelDeltaY && (k = a.wheelDeltaY,
                    c.reverseScrolling && (k = -k), c.tierHolder.scrollTop += k)
            }, !1);
            this.tierHolderHolder.addEventListener("MozMousePixelScroll", function(a) {
                a.stopPropagation();
                a.preventDefault();
                1 == a.axis ? 0 != a.detail && (a = a.detail / 4, c.reverseScrolling && (a = -a), c.move(a)) : (a = a.detail, c.reverseScrolling || (a = -a), c.tierHolder.scrollTop += a)
            }, !1);
            this.tierHolderHolder.addEventListener("touchstart", function(a) {
                return c.touchStartHandler(a)
            }, !1);
            this.tierHolderHolder.addEventListener("touchmove", function(a) {
                return c.touchMoveHandler(a)
            }, !1);
            this.tierHolderHolder.addEventListener("touchend", function(a) {
                return c.touchEndHandler(a)
            }, !1);
            this.tierHolderHolder.addEventListener("touchcancel", function(a) {
                return c.touchCancelHandler(a)
            }, !1);
            var a = function(a) {
                if (13 == a.keyCode) {
                    var k = !1;
                    for (a = 0; a < c.tiers.length; ++a) {
                        var b = c.tiers[a];
                        b.wantedLayoutHeight && b.wantedLayoutHeight != b.layoutHeight && (b.layoutHeight = b.wantedLayoutHeight, b.clipTier(), k = !0)
                    }
                    k && c.arrangeTiers()
                } else if (32 == a.keyCode || 32 == a.charCode) c.isSnapZooming ? (c.isSnapZooming = !1,
                    k = (c.savedZoom || 20) + c.zoomMin) : (c.isSnapZooming = !0, k = (c.savedZoom || 0) + c.zoomMin), c.savedZoom = c.zoomSliderValue - c.zoomMin, c.zoomSliderValue = k, c.zoom(Math.exp(1 * k / c.zoomExpt)), a.stopPropagation(), a.preventDefault();
                else if (85 == a.keyCode) "opts" === c.uiMode && (k = document.getElementById("singleBaseHightlightButton").checked, document.getElementById("singleBaseHightlightButton").checked = !k), c.singleBaseHighlight = !c.singleBaseHighlight, c.positionRuler(), a.stopPropagation(), a.preventDefault();
                else if (39 == a.keyCode) a.stopPropagation(),
                    a.preventDefault(), c.scrollArrowKey(a, -1);
                else if (37 == a.keyCode) a.stopPropagation(), a.preventDefault(), c.scrollArrowKey(a, 1);
                else if (38 == a.keyCode || 87 == a.keyCode)
                    if (a.stopPropagation(), a.preventDefault(), a.shiftKey) a = c.getSelectedTier(), 0 > a || (k = c.tiers[a], a = k.forceHeight || k.subtiers[0].height, 40 <= a && k.mergeConfig({
                        height: a - 10
                    }));
                    else if (a.ctrlKey || a.metaKey) {
                    if (a = c.getSelectedTier(), !(0 > a) && (k = c.tiers[a], k.quantLeapThreshold)) {
                        var b = k.subtiers[0].height,
                            d = k.subtiers[0].quant;
                        d && (a = 1 * d.min, d = 1 * d.max,
                            b = (d - a) / b, k.mergeConfig({
                                quantLeapThreshold: a + ((Math.round((k.quantLeapThreshold - a) / b) | 0) + 1) * b
                            }), k.notify("Threshold: " + Q(k.quantLeapThreshold)))
                    }
                } else if (a.altKey) {
                    if (d = c.selectedTiers.length, 0 != d) {
                        a = c.selectedTiers[0];
                        for (var f = !0, k = [], b = 0; b < c.selectedTiers.length; ++b) k.push(c.tiers[c.selectedTiers[b]]), 0 < b && 1 != c.selectedTiers[b] - c.selectedTiers[b - 1] && (f = !1);
                        if (!(f && 0 >= a)) {
                            for (b = c.selectedTiers.length - 1; 0 <= b; --b) c.tiers.splice(c.selectedTiers[b], 1);
                            c.selectedTiers.splice(0, d);
                            a = f ? a - 1 : a;
                            for (b = 0; b <
                                k.length; ++b) c.tiers.splice(a + b, 0, k[b]), c.selectedTiers.push(a + b);
                            c.withPreservedSelection(c._ensureTiersGrouped);
                            c.markSelectedTiers();
                            c.notifyTierSelection();
                            c.reorderTiers();
                            c.notifyTier()
                        }
                    }
                } else if (a = c.getSelectedTier(), 0 < a) {
                    if (c.setSelectedTier(a - 1), a = c.tiers[c.getSelectedTier()], k = a.row.offsetTop, a = k + a.row.offsetHeight, k < c.tierHolder.scrollTop || a > c.tierHolder.scrollTop + c.tierHolder.offsetHeight) c.tierHolder.scrollTop = k
                } else c.notifyTierSelectionWrap(-1);
                else if (40 == a.keyCode || 83 == a.keyCode)
                    if (a.stopPropagation(),
                        a.preventDefault(), a.shiftKey) a = c.getSelectedTier(), 0 > a || (k = c.tiers[a], a = k.forceHeight || k.subtiers[0].height, k.mergeConfig({
                        height: a + 10
                    }));
                    else if (a.ctrlKey || a.metaKey) {
                    if (a = c.getSelectedTier(), !(0 > a) && (k = c.tiers[a], k.quantLeapThreshold && (b = k.subtiers[0].height, d = k.subtiers[0].quant))) a = 1 * d.min, d = 1 * d.max, b = (d - a) / b, d = Math.round((k.quantLeapThreshold - a) / b) | 0, 1 < d && (k.mergeConfig({
                        quantLeapThreshold: a + (d - 1) * b
                    }), k.notify("Threshold: " + Q(k.quantLeapThreshold)))
                } else if (a.altKey) {
                    if (d = c.selectedTiers.length,
                        0 != d) {
                        a = c.selectedTiers[0];
                        for (var g = 0, k = [], b = 0; b < c.selectedTiers.length; ++b) k.push(c.tiers[c.selectedTiers[b]]), 0 < b && (g += c.selectedTiers[b] - c.selectedTiers[b - 1] - 1);
                        f = 0 == g;
                        if (!(f && a + d >= c.tiers.length)) {
                            for (b = c.selectedTiers.length - 1; 0 <= b; --b) c.tiers.splice(c.selectedTiers[b], 1);
                            c.selectedTiers.splice(0, d);
                            a = f ? a + 1 : a + g;
                            for (b = 0; b < k.length; ++b) c.tiers.splice(a + b, 0, k[b]), c.selectedTiers.push(a + b);
                            c.withPreservedSelection(function() {
                                c._ensureTiersGrouped(!0)
                            });
                            c.markSelectedTiers();
                            c.notifyTierSelection();
                            c.reorderTiers();
                            c.notifyTier()
                        }
                    }
                } else {
                    if (a = c.getSelectedTier(), a < c.tiers.length - 1 && (c.setSelectedTier(a + 1), a = c.tiers[c.getSelectedTier()], k = a.row.offsetTop, a = k + a.row.offsetHeight, k < c.tierHolder.scrollTop || a > c.tierHolder.scrollTop + c.tierHolder.offsetHeight)) c.tierHolder.scrollTop = Math.min(k, a - c.tierHolder.offsetHeight)
                } else if (187 == a.keyCode || 61 == a.keyCode) a.stopPropagation(), a.preventDefault(), c.zoomStep(-10);
                else if (189 == a.keyCode || 173 == a.keyCode) a.stopPropagation(), a.preventDefault(), c.zoomStep(10);
                else if (73 == a.keyCode || 105 == a.keyCode) a.stopPropagation(), a.preventDefault(), a = c.getSelectedTier(), 0 > a || (b = c.tiers[a], b.infoVisible ? (b.infoElement.style.display = "none", b.updateHeight(), b.infoVisible = !1) : (b.infoElement.style.display = "block", b.updateHeight(), b.infoVisible = !0));
                else if (84 == a.keyCode || 116 == a.keyCode)
                    if (a.shiftKey)
                        for (a.stopPropagation(), a.preventDefault(), a = 0; a < c.tiers.length; ++a) b = c.tiers[a], b.dasSource.collapseSuperGroups && (void 0 === k && (k = !b.bumped), b.mergeConfig({
                            bumped: k
                        }));
                    else a.ctrlKey ||
                        a.metaKey || (a.stopPropagation(), a.preventDefault(), a = c.getSelectedTier(), 0 > a || (b = c.tiers[a], b.dasSource.collapseSuperGroups && (void 0 === k && (k = !b.bumped), b.mergeConfig({
                            bumped: k
                        }))));
                else if (77 == a.keyCode || 109 == a.keyCode) a.stopPropagation(), a.preventDefault(), (a.ctrlKey || a.metaKey) && 1 < c.selectedTiers.length && c.mergeSelectedTiers();
                else if (68 == a.keyCode || 100 == a.keyCode) {
                    if (a.stopPropagation(), a.preventDefault(), a.ctrlKey || a.metaKey) a = c.getSelectedTier(), 0 > a || c.addTier(c.tiers[a].dasSource)
                } else if (80 ==
                    a.keyCode || 112 == a.keyCode)
                    if (a.ctrlKey || a.metaKey) {
                        k = [];
                        for (a = 0; a < c.selectedTiers.length; ++a) k.push(c.tiers[c.selectedTiers[a]]);
                        for (a = 0; a < k.length; ++a) k[a].mergeConfig({
                            pinned: !k[a].pinned
                        })
                    }
            };
            this.browserHolder.addEventListener("focus", function(k) {
                c.browserHolder.addEventListener("keydown", a, !1)
            }, !1);
            this.browserHolder.addEventListener("blur", function(k) {
                c.browserHolder.removeEventListener("keydown", a, !1)
            }, !1);
            this.hPopupHolder = f("div");
            this.hPopupHolder.style["font-family"] = "helvetica";
            this.hPopupHolder.style["font-size"] =
                "12pt";
            this.hPopupHolder.classList.add("dalliance");
            document.body.appendChild(this.hPopupHolder);
            for (var b = 0; b < this.sources.length; ++b) {
                var g = this.sources[b];
                if (g) {
                    var k = {};
                    this.restoredConfigs && (k = this.restoredConfigs[b]);
                    g.disabled || this.makeTier(g, k)
                }
            }
            c._ensureTiersGrouped();
            c.arrangeTiers();
            c.reorderTiers();
            c.refresh();
            c.setSelectedTier(1);
            c.positionRuler();
            (b = this.getSequenceSource()) && b.getSeqInfo(this.chr, function(a) {
                c.currentSeqMax = a ? a.length : -1
            });
            this.queryRegistry();
            for (var l in this.chains) this.queryRegistry(l, !0);
            if (this.hubs)
                for (l = 0; l < this.hubs.length; ++l) b = this.hubs[l], "string" == typeof b && (b = {
                        url: b
                    }),
                    function(a) {
                        L(a.url, function(k, b) {
                            if (b) console.log(b);
                            else {
                                var d;
                                if (d = a.genome ? k.genomes[a.genome] : k.genomes[c.coordSystem.ucscName]) a.mapping && (d.mapping = a.mapping), c.hubObjects.push(d)
                            }
                        }, a)
                    }(b);
            this.fullScreen && this.setFullScreenHeight();
            !this.statusRestored && this.storeStatus && this.storeStatus();
            for (l = 0; l < this.initListeners.length; ++l) try {
                this.initListeners[l].call(this)
            } catch (w) {
                console.log(w)
            }
        };
        n.prototype.touchStartHandler =
            function(a) {
                this.touchOriginX = a.touches[0].pageX;
                this.touchOriginY = a.touches[0].pageY;
                2 == a.touches.length && (a = Math.abs(a.touches[0].pageX - a.touches[1].pageX), this.zooming = !0, this.zoomLastSep = this.zoomInitialSep = a, this.zoomInitialScale = this.scale)
            };
        n.prototype.touchMoveHandler = function(a) {
            a.stopPropagation();
            a.preventDefault();
            if (1 == a.touches.length) {
                var c = a.touches[0].pageX;
                a = a.touches[0].pageY;
                this.touchOriginX && c != this.touchOriginX && this.move(c - this.touchOriginX);
                this.touchOriginY && a != this.touchOriginY &&
                    (this.tierHolder.scrollTop -= a - this.touchOriginY);
                this.touchOriginX = c;
                this.touchOriginY = a
            } else if (this.zooming && 2 == a.touches.length) {
                c = Math.abs(a.touches[0].pageX - a.touches[1].pageX);
                if (c != this.zoomLastSep) {
                    a = (a.touches[0].pageX + a.touches[1].pageX) / 2;
                    var b = this.viewStart + a / this.scale | 0;
                    this.scale = c / this.zoomInitialSep * this.zoomInitialScale;
                    this.viewStart = b - a / this.scale | 0;
                    for (a = 0; a < this.tiers.length; ++a) this.tiers[a].draw()
                }
                this.zoomLastSep = c
            }
        };
        n.prototype.touchEndHandler = function(a) {};
        n.prototype.touchCancelHandler =
            function(a) {};
        n.prototype.makeTier = function(a, c) {
            try {
                return this.realMakeTier(a, c)
            } catch (b) {
                console.log("Error initializing", a), console.log(b.stack || b)
            }
        };
        n.prototype.realMakeTier = function(c, b) {
            var d = this,
                f = null;
            this.tierBackgroundColors && (f = this.tierBackgroundColors[this.tiers.length % this.tierBackgroundColors.length]);
            var k = new P(this, c, b, f);
            k.oorigin = this.viewStart;
            var g = !1,
                l, w, q, h = function(c, b) {
                    var f = k.subtiers;
                    if (f) {
                        var g = 0;
                        for (b -= k.padding; g < f.length && b > f[g].height && g < f.length - 1;) b = b - f[g].height -
                            k.padding, ++g;
                        if (!(g >= f.length)) return f = f[g].glyphs, c -= (k.glyphCacheOrigin - d.viewStart) * d.scale, a(f, c, b)
                    }
                },
                t = function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    a = a.clientX;
                    a != w && (d.move(a - w, !0), w = a);
                    d.isDragging = !0
                },
                B = function(a) {
                    window.removeEventListener("mousemove", t, !0);
                    window.removeEventListener("mouseup", B, !0);
                    d.move(a.clientX - w)
                };
            k.viewport.addEventListener("mousedown", function(a) {
                d.browserHolder.focus();
                a.preventDefault();
                k.row.getBoundingClientRect();
                a = a.clientX;
                window.addEventListener("mousemove",
                    t, !0);
                window.addEventListener("mouseup", B, !0);
                l = w = a;
                d.isDragging = !1
            }, !1);
            k.viewport.addEventListener("mousemove", function(a) {
                var c = k.row.getBoundingClientRect(),
                    b = a.clientX - c.left,
                    f = a.clientY - c.top,
                    c = h(b, f);
                k.row.style.cursor = c && 0 < c.length ? "pointer" : "default";
                q && clearTimeout(q);
                g || (q = setTimeout(function() {
                    var c = h(b, f);
                    c && 0 < c.length && d.notifyFeatureHover(a, c[c.length - 1], c, k)
                }, 1E3))
            });
            var r = null;
            k.viewport.addEventListener("mouseup", function(a) {
                var c = k.row.getBoundingClientRect(),
                    b = a.clientX - c.left,
                    c = a.clientY - c.top,
                    f = h(b, c);
                f && 0 < f.length && !d.isDragging && (r ? (clearTimeout(r), r = null, d.featureDoubleClick(f, b, c)) : r = setTimeout(function() {
                    r = null;
                    d.notifyFeature(a, f[f.length - 1], f, k)
                }, 500));
                if (d.isDragging && b != l && k.sequenceSource) {
                    var c = d.viewStart + b / d.scale,
                        g = d.viewStart + l / d.scale;
                    c < g ? (b = c | 0, c = g | 0) : (b = g | 0, c |= 0);
                    d.notifyRegionSelect(d.chr, b, c)
                }
                d.isDragging = !1
            }, !1);
            k.viewport.addEventListener("mouseout", function(a) {
                g = !1
            });
            k.removeButton.addEventListener("click", function(a) {
                a.stopPropagation();
                a.preventDefault();
                for (a = 0; a < d.tiers.length; ++a)
                    if (d.tiers[a] === k) {
                        d.removeTier({
                            index: a
                        });
                        break
                    }
            }, !1);
            k.nameButton.addEventListener("click", function(a) {
                a.stopPropagation();
                a.preventDefault();
                if (a.shiftKey) {
                    a = -1;
                    for (var c = 0; c < d.tiers.length; ++c)
                        if (d.tiers[c] === k) {
                            a = c;
                            break
                        }
                    0 <= a && (c = d.selectedTiers.indexOf(a), 0 <= c ? d.selectedTiers.splice(c, 1) : (d.selectedTiers.push(a), d.selectedTiers.sort()), d.markSelectedTiers(), d.notifyTierSelection(), 0 < d.selectedTiers.length ? d.browserHolder.focus() : d.notifyTierSelectionWrap(-1))
                } else {
                    for (c =
                        0; c < d.tiers.length; ++c)
                        if (d.tiers[c] === k && (d.browserHolder.focus(), 1 != d.selectedTiers.length || d.selectedTiers[0] != c)) {
                            d.setSelectedTier(c);
                            return
                        }
                    k.infoVisible ? (k.infoElement.style.display = "none", k.updateHeight(), k.infoVisible = !1) : (k.infoElement.style.display = "block", k.updateHeight(), k.infoVisible = !0)
                }
            }, !1);
            k.bumpButton.addEventListener("click", function(a) {
                a.stopPropagation();
                a.preventDefault();
                var c;
                k.dasSource.collapseSuperGroups && (void 0 === c && (c = !k.bumped), k.mergeConfig({
                    bumped: c
                }))
            }, !1);
            var e,
                m, D, y, n, H = !1,
                I = function(a) {
                    var c = k.label;
                    a.stopPropagation();
                    a.preventDefault();
                    if (!e) {
                        m = k.pinned ? d.pinnedTierHolder : d.tierHolder;
                        D = m.scrollHeight - m.offsetHeight;
                        e = c.cloneNode(!0);
                        e.style.cursor = "pointer";
                        m.appendChild(e);
                        c.style.visibility = "hidden";
                        for (var b = 0; b < d.tiers.length; ++b)
                            if (d.tiers[b] === k) {
                                y = b;
                                break
                            }
                        n = a.clientY
                    }
                    var f = m.getBoundingClientRect();
                    e.style.left = c.getBoundingClientRect().left - f.left + "px";
                    e.style.top = a.clientY - f.top + m.scrollTop - 10 + "px";
                    c = a.clientY - f.top + m.scrollTop;
                    for (b = 0; b < d.tiers.length; ++b)
                        if (f =
                            d.tiers[b], !(f.pinned ^ k.pinned) && (f = f.row.getBoundingClientRect(), c -= f.bottom - f.top, 0 > c)) {
                            if (b < y && a.clientY < n || b > y && a.clientY > n) {
                                d.withPreservedSelection(function() {
                                    d.tiers.splice(y, 1);
                                    d.tiers.splice(b, 0, k);
                                    d._ensureTiersGrouped(b > y)
                                });
                                for (c = 0; c < d.tiers.length; ++c) d.tiers[c] == k && (y = c);
                                n = a.clientY;
                                d.reorderTiers();
                                m.appendChild(e);
                                H = !0
                            }
                            break
                        }
                    e.offsetTop < m.scrollTop ? m.scrollTop -= m.scrollTop - e.offsetTop : e.offsetTop + e.offsetHeight > m.scrollTop + m.offsetHeight && (m.scrollTop = Math.min(m.scrollTop + (e.offsetTop +
                        e.offsetHeight) - (m.scrollTop + m.offsetHeight), D))
                },
                p = function(a) {
                    var c = k.label;
                    a.stopPropagation();
                    a.preventDefault();
                    e && (e.style.cursor = "auto", m.removeChild(e), e = null, c.style.visibility = "visible");
                    document.removeEventListener("mousemove", I, !1);
                    document.removeEventListener("mouseup", p, !1);
                    if (H) {
                        for (a = 0; a < d.tiers.length; ++a)
                            if (d.tiers[a] == k) {
                                d.setSelectedTier(a);
                                break
                            }
                        d.notifyTier()
                    }
                };
            k.label.addEventListener("mousedown", function(a) {
                a.stopPropagation();
                a.preventDefault();
                H = !1;
                document.addEventListener("mousemove",
                    I, !1);
                document.addEventListener("mouseup", p, !1)
            }, !1);
            this.tiers.push(k);
            k.init();
            k.currentlyHeight = 50;
            this.updateHeight();
            k.updateLabel();
            k.featureSource && k.featureSource.addActivityListener && k.featureSource.addActivityListener(function(a) {
                k.loaderButton.style.display = 0 < a ? "inline-block" : "none";
                d.pingActivity()
            });
            this.withPreservedSelection(d._ensureTiersGrouped);
            k._updateFromConfig();
            this.reorderTiers();
            return k
        };
        n.prototype.reorderTiers = function() {
            d(this.tierHolder);
            d(this.pinnedTierHolder);
            this.disablePinning &&
                (this.tierHolder.appendChild(this.ruler), this.tierHolder.appendChild(this.ruler2));
            for (var a = !1, c = [], b = [], f = 0; f < this.tiers.length; ++f) {
                var k = this.tiers[f];
                k.pinned && !this.disablePinning ? (c.push(k), this.pinnedTierHolder.appendChild(this.tiers[f].row), a = !0) : (b.push(k), this.tierHolder.appendChild(this.tiers[f].row))
            }
            this.withPreservedSelection(function() {
                this.tiers.splice(0, this.tiers.length);
                for (var a = 0; a < c.length; ++a) this.tiers.push(c[a]);
                for (a = 0; a < b.length; ++a) this.tiers.push(b[a])
            });
            a ? this.pinnedTierHolder.classList.add("tier-holder-pinned-full") :
                this.pinnedTierHolder.classList.remove("tier-holder-pinned-full");
            this.arrangeTiers()
        };
        n.prototype.withPreservedSelection = function(a) {
            for (var c = [], b = 0; b < this.selectedTiers.length; ++b) c.push(this.tiers[this.selectedTiers[b]]);
            a.call(this);
            this.selectedTiers = [];
            for (a = 0; a < this.tiers.length; ++a) 0 <= c.indexOf(this.tiers[a]) && this.selectedTiers.push(a)
        };
        n.prototype.refreshTier = function(a) {
            this.knownSpace && this.knownSpace.invalidate(a)
        };
        n.prototype._ensureTiersGrouped = function(a) {
            for (var c = {}, b = 0; b < this.tiers.length; ++b) {
                var d =
                    this.tiers[b];
                d.dasSource.tierGroup && pusho(c, d.dasSource.tierGroup, d)
            }
            var k = [];
            a && this.tiers.reverse();
            for (b = 0; b < this.tiers.length; ++b)
                if (d = this.tiers[b], d.dasSource.tierGroup) {
                    var f = c[d.dasSource.tierGroup];
                    if (f) {
                        a && f.reverse();
                        for (var g = 0; g < f.length; ++g) k.push(f[g]);
                        c[d.dasSource.tierGroup] = null
                    }
                } else k.push(d);
            a && k.reverse();
            this.tiers.splice(0, this.tiers.length);
            for (g = 0; g < k.length; ++g) this.tiers.push(k[g])
        };
        n.prototype.arrangeTiers = function() {
            for (var a = [], c = {}, b = 0; b < this.tiers.length; ++b) {
                var d =
                    this.tiers[b];
                d.pinned && (a.push(d), d.dasSource.tierGroup && pusho(c, d.dasSource.tierGroup, d))
            }
            for (b = 0; b < this.tiers.length; ++b) d = this.tiers[b], d.pinned || (a.push(d), d.dasSource.tierGroup && pusho(c, d.dasSource.tierGroup, d));
            for (var k in c) {
                var d = c[k],
                    g = this.tierGroups[k];
                g || (g = {
                    element: f("div", f("span", k, {
                        className: "tier-group-label"
                    }), {
                        className: "tier-group"
                    })
                }, this.tierGroups[k] = g);
                g.element.parentNode && g.element.parentNode.removeChild(g.element);
                for (var l = d[0].pinned ? this.pinnedTierHolder : this.tierHolder,
                        w = 1E7, q = 0, b = 0; b < d.length; ++b) var h = d[b].row,
                    w = Math.min(w, h.offsetTop),
                    q = Math.max(q, h.offsetTop + h.offsetHeight);
                g.element.style.top = w + "px";
                g.element.style.left = "0px";
                g.element.style.height = q - w + "px";
                l.appendChild(g.element)
            }
            if (this.tierBackgroundColors)
                for (b = 0; b < a.length; ++b) d = a[b], d.setBackground(this.tierBackgroundColors[b % this.tierBackgroundColors.length]), d.label.style.left = d.dasSource.tierGroup ? "18px" : "2px", d.background = this.tierBackgroundColors[b % this.tierBackgroundColors.length]
        };
        n.prototype.refresh =
            function() {
                this.notifyLocation();
                var a = 100 / this.scale | 0,
                    c = 1E3 / this.scale | 0,
                    b = (this.viewStart + this.viewEnd) / 2,
                    d = b - this.origin;
                this.origin = b;
                this.scaleAtLastRedraw = this.scale;
                for (b = 0; b < this.tiers.length; ++b) {
                    var k = d;
                    this.tiers[b].originHaxx && (k += this.tiers[b].originHaxx);
                    this.tiers[b].originHaxx = k
                }
                d = this.targetQuantRes / this.scale;
                b = Math.max(1, (this.viewStart | 0) - a);
                a = Math.min((this.viewEnd | 0) + a, 0 < (this.currentSeqMax | 0) ? this.currentSeqMax | 0 : 1E9);
                k = Math.max(1, (this.viewStart | 0) - c);
                c = Math.min((this.viewEnd |
                    0) + c, 0 < (this.currentSeqMax | 0) ? this.currentSeqMax | 0 : 1E9);
                if (!this.knownSpace || this.knownSpace.chr !== this.chr) {
                    var f = this.getSequenceSource();
                    this.knownSpace && this.knownSpace.cancel();
                    this.knownSpace = new B(this.tiers, this.chr, k, c, d, f)
                }(f = this.knownSpace.bestCacheOverlapping(this.chr, b, a)) && f.min <= b && f.max >= a ? (this.drawnStart = Math.max(f.min, k), this.drawnEnd = Math.min(f.max, c)) : (this.drawnStart = k, this.drawnEnd = c);
                this.knownSpace.viewFeatures(this.chr, this.drawnStart, this.drawnEnd, d);
                this.drawOverlays();
                this.positionRuler()
            };
        n.prototype.queryRegistry = function(a, c) {
            var b, d;
            a ? (b = this.chains[a].coords, this.mappableSources[a] || (this.mappableSources[a] = new l), d = this.mappableSources[a]) : (b = this.coordSystem, d = this.availableSources);
            var k = G(t(b));
            if (c) {
                var f = localStorage["dalliance.registry." + k + ".last_queried"];
                if (f) try {
                    if (m(d, JSON.parse(localStorage["dalliance.registry." + k + ".sources"]), a), 432E5 > (Date.now() | 0) - (f | 0)) return
                } catch (g) {
                    console.log("Bad registry cache: " + g)
                }
            }
            f = this.registry;
            if (0 == f.indexOf("//")) {
                var w =
                    window.location.protocol;
                "https:" != w && "http:" != w && (f = "http:" + f)
            }(new I(f)).sources(function(c) {
                for (var f = [], g = 0; g < c.length; ++g) {
                    var l = c[g];
                    if (l.coords && 0 != l.coords.length) {
                        var w = l.coords[0];
                        w.taxon == b.taxon && w.auth == b.auth && w.version == b.version && f.push(l)
                    }
                }
                localStorage["dalliance.registry." + k + ".sources"] = JSON.stringify(f);
                localStorage["dalliance.registry." + k + ".last_queried"] = "" + Date.now();
                m(d, f, a)
            }, function(a) {}, b)
        };
        n.prototype.move = function(a, c) {
            var b = this.viewEnd - this.viewStart,
                d = this.viewStart -
                1 * a / this.scale,
                k = d + b;
            c || (0 < this.currentSeqMax && k > this.currentSeqMax && (k = this.currentSeqMax, d = this.viewEnd - b), 1 > d && (d = 1, k = d + b));
            this.setLocation(null, d, k, null, c)
        };
        n.prototype.zoomStep = function(a) {
            var c = 1 * this.zoomSliderValue;
            a = c + a;
            a < this.zoomMin && (a = this.zoomMin);
            a > this.zoomMax && (a = this.zoomMax);
            a != c && (this.zoomSliderValue = a, this.zoom(Math.exp(1 * a / this.zoomExpt)))
        };
        n.prototype.zoom = function(a) {
            this.zoomFactor = a;
            a = Math.round((this.viewStart + this.viewEnd) / 2) | 0;
            this.viewStart = a - this.zoomBase * this.zoomFactor /
                2;
            this.viewEnd = a + this.zoomBase * this.zoomFactor / 2;
            0 < this.currentSeqMax && this.viewEnd > this.currentSeqMax + 5 && (a = this.viewEnd - this.viewStart + 1, this.viewEnd = this.currentSeqMax, this.viewStart = this.viewEnd - a + 1);
            1 > this.viewStart && (a = this.viewEnd - this.viewStart + 1, this.viewStart = 1, this.viewEnd = this.viewStart + a - 1);
            this.scale = this.featurePanelWidth / (this.viewEnd - this.viewStart);
            this.notifyLocation();
            this.refresh()
        };
        n.prototype.spaceCheck = function(a) {
            this.knownSpace && this.knownSpace.chr === this.chr ? (a = 100 / this.scale |
                0, ((this.drawnStart | 0) > Math.max(1, (this.viewStart | 0) - a | 0) || (this.drawnEnd | 0) < Math.min((this.viewEnd | 0) + a, 0 < (this.currentSeqMax | 0) ? this.currentSeqMax | 0 : 1E9)) && this.refresh()) : this.refresh()
        };
        n.prototype.resizeViewer = function(a) {
            var c = this.tierHolder.getBoundingClientRect().width | 0;
            if (0 != c) {
                var b = Math.max(this.featurePanelWidth, 300);
                this.featurePanelWidth = c | 0;
                b != this.featurePanelWidth && (this.zoomMax = this.zoomExpt * Math.log(this.maxViewWidth / this.zoomBase), this.zoomMin = this.zoomExpt * Math.log(this.featurePanelWidth /
                    this.maxPixelsPerBase / this.zoomBase), this.zoomSliderValue = this.zoomExpt * Math.log((this.viewEnd - this.viewStart + 1) / this.zoomBase), this.viewEnd = this.viewStart + (this.viewEnd - this.viewStart) * this.featurePanelWidth / b, c = this.viewEnd - this.viewStart + 1, 0 < this.currentSeqMax && this.viewEnd > this.currentSeqMax && (this.viewEnd = this.currentSeqMax, this.viewStart = this.viewEnd - c + 1), 1 > this.viewStart && (this.viewStart = 1, this.viewEnd = this.viewStart + c - 1), this.positionRuler(), a || this.spaceCheck(), this.notifyLocation());
                this.fullScreen &&
                    this.setFullScreenHeight()
            }
        };
        n.prototype.setFullScreenHeight = function() {
            this.browserHolder.style.maxHeight = Math.max(300, window.innerHeight - (document.body.offsetHeight - this.browserHolder.offsetHeight) - 20) + "px"
        };
        n.prototype.addTier = function(a) {
            a = D(a);
            a.disabled = !1;
            a = this.makeTier(a);
            this.markSelectedTiers();
            this.positionRuler();
            this.notifyTier();
            return a
        };
        n.prototype.removeTier = function(a, c) {
            var b = -1;
            if ("undefined" !== typeof a.index && 0 <= a.index && a.index < this.tiers.length) b = a.index;
            else
                for (var d = 0; d <
                    this.tiers.length; ++d)
                    if (z(a, this.tiers[d].dasSource)) {
                        b = d;
                        break
                    } if (0 > b) throw "Couldn't find requested tier";
            this.tiers.splice(b, 1);
            for (var d = [], k = 0; k < this.selectedTiers.length; ++k) {
                var f = this.selectedTiers[k];
                f < b ? d.push(f) : f > b && d.push(f - 1)
            }
            this.selectedTiers = d;
            this.markSelectedTiers();
            this.reorderTiers();
            this.notifyTier()
        };
        n.prototype.getSequenceSource = function() {
            void 0 === this._sequenceSource && (this._sequenceSource = this._getSequenceSource());
            return this._sequenceSource
        };
        n.prototype._getSequenceSource =
            function() {
                for (var a = 0; a < this.tiers.length; ++a)
                    if (this.tiers[a].sequenceSource) return this.tiers[a].sequenceSource;
                for (a = 0; a < this.defaultSources.length; ++a) {
                    var b = this.defaultSources[a];
                    if (b.provides_entrypoints || "sequence" == b.tier_type || b.twoBitURI || b.twoBitBlob) return b.twoBitURI || b.twoBitBlob ? new c(b) : new w(b)
                }
            };
        n.prototype.setLocation = function(a, c, b, d, k) {
            if ("number" !== typeof c) throw Error("minimum must be a number (got " + JSON.stringify(c) + ")");
            if ("number" !== typeof b) throw Error("maximum must be a number (got " +
                JSON.stringify(b) + ")");
            if (c > b) {
                var f = c;
                c = b;
                b = f
            } else c === b && (b += 1);
            d || (d = function(a) {
                if (a) throw a;
            });
            var g = this;
            if ((!a || a == this.chr) && 0 < this.currentSeqMax) return this._setLocation(null, c, b, null, d, k);
            var l = this.getSequenceSource();
            if (!l) return d("Need a sequence source");
            var w = a || this.chr;
            l.getSeqInfo(w, function(f) {
                if (f) return g._setLocation(a, c, b, f, d, k);
                var q;
                q = 0 == w.indexOf("chr") ? w.substr(3) : "chr" + w;
                l.getSeqInfo(q, function(f) {
                    return !f && a ? d("Couldn't find sequence '" + a + "'") : f ? g._setLocation(q, c, b,
                        f, d, k) : g._setLocation(null, c, b, null, d, k)
                })
            })
        };
        n.prototype._setLocation = function(a, c, b, d, k, f) {
            var g = !1;
            a && (0 == a.indexOf("chr") && (a = a.substring(3)), this.chr != a && (g = !0), this.chr = a, this.currentSeqMax = d.length);
            c = parseFloat(c);
            b = parseFloat(b);
            a = Math.max(10, b - c + 1);
            f || (f = this.currentSeqMax, 0 >= f && (f = 1E12), 1 > c && (c = 1, b = c + a - 1), b > f && (b = f, c = Math.max(1, b - a + 1)));
            this.viewStart = c;
            this.viewEnd = b;
            c = Math.max(this.featurePanelWidth, 50) / (this.viewEnd - this.viewStart);
            f = 1E-6 < Math.abs(c - this.scale);
            this.scale = c;
            b = this.zoomSliderValue;
            this.zoomSliderValue = c = this.zoomExpt * Math.log((this.viewEnd - this.viewStart + 1) / this.zoomBase);
            if (f || g) {
                for (g = 0; g < this.tiers.length; ++g) this.tiers[g].viewportHolder.style.left = "5000px", this.tiers[g].overlay.style.left = "5000px";
                this.refresh();
                this.savedZoom ? (c -= this.zoomMin, b -= this.zoomMin, g = c - this.savedZoom, Math.abs(c - b) > Math.abs(g) && (this.isSnapZooming = !this.isSnapZooming, this.savedZoom = b)) : (this.isSnapZooming = !1, this.savedZoom = null)
            } else
                for (g = 0; g < this.tiers.length; ++g) this.tiers[g].viewportHolder.style.left =
                    "" + ((-((this.viewStart - this.tiers[g].norigin) * this.scale) | 0) - 1E3) + "px", this.tiers[g].drawOverlay();
            this.notifyLocation();
            this.spaceCheck();
            this.instrumentActivity && (this.activityStartTime = Date.now() | 0);
            return k()
        };
        n.prototype.setCenterLocation = function(a, c) {
            var b = (this.viewEnd - this.viewStart) / 2;
            this.setLocation(a, c - b, c + b)
        };
        n.prototype.pingActivity = function() {
            if (this.instrumentActivity && this.activityStartTime) {
                for (var a = 0, c = 0; c < this.tiers.length; ++c) "none" !== this.tiers[c].loaderButton.style.display &&
                    ++a;
                0 == a && (a = Date.now() | 0, console.log("Loading took " + (a - this.activityStartTime) + "ms"), this.activityStartTime = null)
            }
        };
        n.prototype.addInitListener = function(a) {
            this.initListeners.push(a)
        };
        n.prototype.addFeatureListener = function(a, c) {
            this.featureListeners.push(a)
        };
        n.prototype.notifyFeature = function(a, c, b, d) {
            for (var k = 0; k < this.featureListeners.length; ++k) try {
                if (this.featureListeners[k](a, c, b, d)) break
            } catch (f) {
                console.log(f.stack)
            }
        };
        n.prototype.addFeatureHoverListener = function(a, c) {
            this.featureHoverListeners.push(a)
        };
        n.prototype.notifyFeatureHover = function(a, c, b, d) {
            for (var k = 0; k < this.featureHoverListeners.length; ++k) try {
                this.featureHoverListeners[k](a, c, b, d)
            } catch (f) {
                console.log(f.stack)
            }
        };
        n.prototype.addViewListener = function(a, c) {
            this.viewListeners.push(a)
        };
        n.prototype.notifyLocation = function() {
            var a = Math.max(1, this.viewStart | 0),
                c = this.viewEnd | 0;
            0 < this.currentSeqMax && c > this.currentSeqMax && (c = this.currentSeqMax);
            for (var b = 0; b < this.viewListeners.length; ++b) try {
                this.viewListeners[b](this.chr, a, c, this.zoomSliderValue, {
                    current: this.zoomSliderValue,
                    alternate: this.savedZoom + this.zoomMin || this.zoomMin,
                    isSnapZooming: this.isSnapZooming,
                    min: this.zoomMin,
                    max: this.zoomMax
                }, this.viewStart, this.viewEnd)
            } catch (d) {
                console.log(d.stack)
            }
        };
        n.prototype.addTierListener = function(a) {
            this.tierListeners.push(a)
        };
        n.prototype.notifyTier = function() {
            for (var a = 0; a < this.tierListeners.length; ++a) try {
                this.tierListeners[a]()
            } catch (c) {
                console.log(c.stack)
            }
        };
        n.prototype.addRegionSelectListener = function(a) {
            this.regionSelectListeners.push(a)
        };
        n.prototype.notifyRegionSelect = function(a, c, b) {
            for (var d = 0; d < this.regionSelectListeners.length; ++d) try {
                this.regionSelectListeners[d](a, c, b)
            } catch (k) {
                console.log(k.stack)
            }
        };
        n.prototype.highlightRegion = function(a, c, b) {
            var d = this;
            if (a == this.chr) return this._highlightRegion(a, c, b);
            var k = this.getSequenceSource();
            if (!k) throw "Need a sequence source";
            k.getSeqInfo(a, function(f) {
                if (f) return d._highlightRegion(a, c, b);
                var g;
                g = 0 == a.indexOf("chr") ? a.substr(3) : "chr" + a;
                k.getSeqInfo(g, function(a) {
                    if (a) return d._highlightRegion(g,
                        c, b)
                })
            })
        };
        n.prototype._highlightRegion = function(a, c, b) {
            for (var d = 0; d < this.highlights.length; ++d) {
                var k = this.highlights[d];
                if (k.chr == a && k.min == c && k.max == b) return
            }
            this.highlights.push(new p(a, c, b));
            d = this.viewStart - 1E3 / this.scale;
            k = this.viewEnd + 1E3 / this.scale;
            (a == this.chr || a == "chr" + this.chr) && c < k && b > d && this.drawOverlays();
            this.notifyLocation()
        };
        n.prototype.clearHighlights = function() {
            this.highlights = [];
            this.drawOverlays();
            this.notifyLocation()
        };
        n.prototype.drawOverlays = function() {
            for (var a = 0; a < this.tiers.length; ++a) this.tiers[a].drawOverlay()
        };
        n.prototype.featuresInRegion = function(a, c, b) {
            var d = [];
            if (a !== this.chr) return [];
            for (a = 0; a < this.tiers.length; ++a)
                for (var k = this.tiers[a].currentFeatures || [], f = 0; f < k.length; ++f) {
                    var g = k[f];
                    g.min <= b && g.max >= c && d.push(g)
                }
            return d
        };
        n.prototype.getSelectedTier = function() {
            return 0 < this.selectedTiers.length ? this.selectedTiers[0] : -1
        };
        n.prototype.setSelectedTier = function(a) {
            this.selectedTiers = null == a ? [] : [a];
            this.markSelectedTiers();
            this.notifyTierSelection()
        };
        n.prototype.markSelectedTiers = function() {
            for (var a =
                    0; a < this.tiers.length; ++a) {
                var c = this.tiers[a].nameButton;
                0 <= this.selectedTiers.indexOf(a) ? c.classList.add("active") : c.classList.remove("active")
            }
            0 < this.selectedTiers.length && (a = this.browserHolder.offsetTop + this.browserHolder.offsetHeight / 2, a > document.body.scrollTop && a + 100 < document.body.scrollTop + window.innerHeight && this.browserHolder.focus())
        };
        n.prototype.addTierSelectionListener = function(a) {
            this.tierSelectionListeners.push(a)
        };
        n.prototype.notifyTierSelection = function() {
            for (var a = 0; a < this.tierSelectionListeners.length; ++a) try {
                this.tierSelectionListeners[a](this.selectedTiers)
            } catch (c) {
                console.log(c.stack)
            }
        };
        n.prototype.addTierSelectionWrapListener = function(a) {
            this.tierSelectionWrapListeners.push(a)
        };
        n.prototype.notifyTierSelectionWrap = function(a) {
            for (var c = 0; c < this.tierSelectionWrapListeners.length; ++c) try {
                this.tierSelectionWrapListeners[c](a)
            } catch (b) {
                console.log(b.stack)
            }
        };
        n.prototype.positionRuler = function() {
            var a = "none",
                c = "",
                b = "";
            "center" == this.rulerLocation ? (a = "block", c = "" + (this.featurePanelWidth / 2 | 0) + "px") : "left" == this.rulerLocation ? (a = "block", c = "0px") : "right" == this.rulerLocation ? (a = "block", b =
                "0px") : a = "none";
            this.ruler.style.display = a;
            this.ruler.style.left = c;
            this.ruler.style.right = b;
            this.singleBaseHighlight ? (this.ruler2.style.display = "block", this.ruler2.style.borderWidth = "1px", 1 > this.scale ? (this.ruler2.style.width = "0px", this.ruler2.style.borderRightWidth = "0px") : (this.ruler2.style.width = this.scale + "px", this.ruler2.style.borderRightWidth = "1px"), this.locSingleBase.style.visibility = "visible", this.locSingleBase.style.left = "" + (this.featurePanelWidth / 2 - this.locSingleBase.offsetWidth / 2 + this.ruler2.offsetWidth /
                2 | 0) + "px") : (this.locSingleBase.style.visibility = "hidden", this.ruler2.style.width = "1px", this.ruler2.style.borderWidth = "0px", this.ruler2.style.display = "center" == this.rulerLocation ? "none" : "block");
            this.ruler2.style.left = "" + (this.featurePanelWidth / 2 | 0) + "px";
            for (var d = 0; d < this.tiers.length; ++d) {
                var k = this.tiers[d],
                    f = k.quantOverlay,
                    g;
                k.subtiers && 0 < k.subtiers.length && (g = k.subtiers[0].quant);
                f && (f.style.display = g ? a : "none", f.style.left = c, f.style.right = b)
            }
        };
        n.prototype.featureDoubleClick = function(a, c, b) {
            if (a &&
                0 != a.length && (b = a[a.length - 1], b.min && b.max)) {
                var d = ((b.min | 0) - (this.viewStart | 0)) * this.scale,
                    k = (b.max - b.min + 1) * this.scale;
                a = ((b.min | 0) + (b.max | 0)) / 2;
                10 < k && (c = 1 * (c - d) / k, 0.3 > c ? a = b.min | 0 : 0.7 < c && (a = (b.max | 0) + 1));
                c = this.viewEnd - this.viewStart;
                this.setLocation(null, a - c / 2, a + c / 2)
            }
        };
        n.prototype.zoomForScale = function(a) {
            return 0.2 < a ? "high" : 0.01 < a ? "medium" : "low"
        };
        n.prototype.zoomForCurrentScale = function() {
            return this.zoomForScale(this.scale)
        };
        n.prototype.updateHeight = function() {
            for (var a = 0, c = 0; c < this.tiers.length; ++c) a +=
                this.tiers[c].currentHeight || 30;
            this.ruler.style.height = "" + a + "px";
            this.ruler2.style.height = "" + a + "px";
            this.browserHolder.style.display = "block";
            this.browserHolder.style.display = "flex"
        };
        n.prototype.scrollArrowKey = function(a, c) {
            this.reverseKeyScrolling && (c = -c);
            if (a.ctrlKey || a.metaKey) {
                var b = !1;
                a.shiftKey && (b = !0);
                this.leap(c, b)
            } else if (1 < this.scale) {
                var b = (this.viewStart + this.viewEnd) / 2,
                    b = b - Math.round(b),
                    d = 1;
                a.shiftKey && (d *= 10);
                0 < c ? (d = -d - b, 0 < b && (d += 1)) : (d -= b, 0 > b && (d -= 1));
                this.setLocation(null, this.viewStart +
                    d, this.viewEnd + d)
            } else this.move(a.shiftKey ? 100 * c : 25 * c)
        };
        n.prototype.leap = function(a, c) {
            var b = this,
                d = (b.viewStart + b.viewEnd + 1) / 2 | 0;
            0 < a && 1 >= b.viewStart ? d -= 1E8 : 0 > a && b.viewEnd >= b.currentSeqMax && (d += 1E8);
            var k = b.getSelectedTier();
            0 > k || ((k = b.tiers[k]) && (k.featureSource && this.sourceAdapterIsCapable(k.featureSource, "quantLeap") && "number" == typeof k.quantLeapThreshold || k.featureSource && this.sourceAdapterIsCapable(k.featureSource, "leap")) ? k.findNextFeature(b.chr, d, -a, c, function(k) {
                if (k) {
                    var f = k.min,
                        g = k.max;
                    c && (0 < a ? f > d + 1 ? g = f : (g++, f = g) : g < d - 1 ? (g++, f = g) : g = f);
                    var l = b.viewEnd - b.viewStart + 1;
                    parseFloat(l / 2) == parseInt(l / 2) && l--;
                    f = (f + g - l) / 2 + 1;
                    b.setLocation(k.segment, f, f + l - 1)
                } else alert("no next feature")
            }) : this.move(100 * a))
        };
        n.prototype.nameForCoordSystem = function(a) {
            var c = null,
                b = null;
            this.assemblyNamePrimary && (c = "" + a.auth, "undefined" !== typeof a.version && (c += a.version));
            this.assemblyNameUcsc && (b = a.ucscName);
            return null != c && null != b ? c + "/" + b : c || b || "unknown"
        };
        n.prototype.makeLoader = function(a) {
            var c = 1 < window.devicePixelRatio;
            return 20 > (a || 16) ? f("img", null, {
                src: this.resolveURL("$$img/spinner_" + (c ? 16 : 32) + ".gif"),
                width: "16",
                height: "16"
            }) : f("img", null, {
                src: this.resolveURL("$$img/spinner_" + (c ? 24 : 48) + ".gif"),
                width: "24",
                height: "24"
            })
        };
        n.prototype.getWorker = function() {
            if (!this.useFetchWorkers || !this.fetchWorkers || 0 == this.fetchWorkers.length) return null;
            this.nextWorker >= this.fetchWorkers.length && (this.nextWorker = 0);
            return this.fetchWorkers[this.nextWorker++]
        };
        b.prototype.postCommand = function(a, c, b) {
            var d = "x" + ++this.tagSeed;
            a.tag =
                d;
            this.callbacks[d] = c;
            this.worker.postMessage(a, b)
        };
        if ("undefined" !== typeof u) {
            u.exports = {
                Browser: n,
                sourcesAreEqual: z,
                sourceDataURI: h
            };
            e("./browser-ui");
            e("./track-adder");
            e("./feature-popup");
            e("./tier-actions");
            e("./domui");
            e("./search");
            u = e("./sourceadapters");
            var c = u.TwoBitSequenceSource,
                w = u.DASSequenceSource,
                B = e("./kspace").KnownSpace,
                I = e("./das").DASRegistry
        }
        g.prototype.get = function(a) {
            var c = this.sourcesByURI[h(a)];
            if (c)
                for (var b = 0; b < c.configs.length; ++b)
                    if (r(c.configs[b], a)) return c.sources[b]
        };
        g.prototype.put = function(a, c) {
            var b = h(a),
                d = this.sourcesByURI[b];
            d || (d = {
                configs: [],
                sources: []
            }, this.sourcesByURI[b] = d);
            d.configs.push(a);
            d.sources.push(c)
        }
    }, {
        "./browser-ui": 5,
        "./chainset": 7,
        "./das": 10,
        "./domui": 11,
        "./feature-popup": 19,
        "./kspace": 23,
        "./numformats": 26,
        "./search": 30,
        "./sha1": 33,
        "./sourceadapters": 34,
        "./thub": 41,
        "./tier": 44,
        "./tier-actions": 42,
        "./track-adder": 45,
        "./utils": 48,
        "./version": 50,
        "es6-promise": 52
    }],
    7: [function(e, u, s) {
        function p(a, b, d, f) {
            "string" == typeof a ? (this.uri = a, this.srcTag =
                b, this.destTag = d, this.coords = f) : (this.uri = a.uri, this.srcTag = a.srcTag, this.destTag = a.destTag, this.coords = g(a.coords), this.type = a.type, this.credentials = a.credentials);
            this.chainsBySrc = {};
            this.chainsByDest = {};
            this.postFetchQueues = {};
            this.fetchedTiles = {};
            this.granularity = 1E6;
            this.chainFetcher = "bigbed" == this.type ? new m(this.uri, this.credentials) : "alias" == this.type ? new z(a) : new n(this.uri, this.srcTag, this.destTag)
        }

        function n(b, d, f) {
            this.source = new a(b);
            this.srcTag = d;
            this.destTag = f
        }

        function m(a, b) {
            var g =
                this;
            this.uri = a;
            this.credentials = b;
            this.bwg = new t(function(a, b) {
                d(new f(g.uri, {
                    credentials: g.credentials
                }), function(d, f) {
                    d ? a(d) : b(f)
                })
            });
            this.bwg.then(function(a, b) {
                b && console.log(b)
            })
        }

        function h(a) {
            return parseInt(a)
        }

        function v(a) {
            return 0 == a.indexOf("chr") ? a.substr(3) : a
        }

        function r(a) {
            var b = {
                    srcChr: v(a.srcChrom),
                    srcMin: parseInt(a.srcStart),
                    srcMax: parseInt(a.srcEnd),
                    srcOri: a.srcOri,
                    destChr: v(a.segment),
                    destMin: a.min - 1,
                    destMax: a.max,
                    destOri: a.ori,
                    blocks: []
                },
                d = a.srcStarts.split(",").map(h),
                f = a.destStarts.split(",").map(h);
            a = a.blockLens.split(",").map(h);
            for (var g = 0; g < d.length; ++g) b.blocks.push([d[g], f[g], a[g]]);
            return b
        }

        function z(a) {
            this.conf = a;
            this.forwardAliases = {};
            a = a.sequenceAliases || [];
            for (var b = 0; b < a.length; ++b) {
                var d = a[b];
                if (!(2 > d.length)) {
                    for (var f = [], g = 0; g < d.length - 1; ++g) f.push(d[g]);
                    this.forwardAliases[d[d.length - 1]] = f
                }
            }
        }
        if ("undefined" !== typeof e) {
            s = e("./das");
            var a = s.DASSource,
                b = s.DASSegment;
            s = e("./utils");
            var q = s.pusho,
                g = s.shallowCopy,
                l = e("./cigar").parseCigar,
                f = e("./bin").URLFetchable,
                d = e("./bigwig").makeBwg,
                t = e("es6-promise").Promise
        }
        p.prototype.exportConfig = function() {
            return {
                uri: this.uri,
                srcTag: this.srcTag,
                destTag: this.destTag,
                coords: this.coords,
                type: this.type,
                credentials: this.credentials
            }
        };
        p.prototype.mapPoint = function(a, b) {
            for (var d = this.chainsBySrc[a] || [], f = 0; f < d.length; ++f) {
                var g = d[f];
                if (b >= g.srcMin && b <= g.srcMax) {
                    var l;
                    l = "-" == g.srcOri ? g.srcMax - b : b - g.srcMin;
                    for (var q = g.blocks, h = 0; h < q.length; ++h) {
                        var t = q[h],
                            c = t[0],
                            w = t[1],
                            t = t[2];
                        if (l >= c && l <= c + t) return d = l - c, {
                            seq: g.destChr,
                            pos: "-" == g.destOri ? g.destMax -
                                w - d : d + w + g.destMin,
                            flipped: g.srcOri != g.destOri
                        }
                    }
                }
            }
            return null
        };
        p.prototype.mapSegment = function(a, b, d) {
            a = this.chainsBySrc[a] || [];
            for (var f = [], g = 0; g < a.length; ++g) {
                var l = a[g];
                if (d >= l.srcMin && b <= l.srcMax) {
                    var q, h;
                    "-" == l.srcOri ? (q = l.srcMax - d, h = l.srcMax - b) : (q = b - l.srcMin, h = d - l.srcMin);
                    for (var t = l.blocks, c = 0; c < t.length; ++c) {
                        var w = t[c],
                            B = w[0],
                            e = w[1],
                            w = w[2];
                        if (h >= B && q <= B + w) {
                            var x = {
                                segment: l.destChr,
                                flipped: "-" == l.srcOri ^ "-" == l.destOri
                            };
                            "-" == l.destOri ? (q >= B ? x.max = l.destMax - e - q + B : (x.max = l.destMax - e, x.partialMax =
                                B - q), h <= B + w ? x.min = l.destMax - e - h + B : (x.min = l.destMax - e - w, x.partialMin = h - B - w)) : (q >= B ? x.min = l.destMin + e + q - B : (x.min = l.destMin + e, x.partialMin = B - q), h <= B + w ? x.max = l.destMin + e + h - B : (x.max = l.destMin + e + w, x.partialMax = h - B - w));
                            f.push(x)
                        }
                    }
                }
            }
            return f
        };
        p.prototype.unmapPoint = function(a, b) {
            for (var d = this.chainsByDest[a] || [], f = 0; f < d.length; ++f) {
                var g = d[f];
                if (b >= g.destMin && b <= g.destMax) {
                    var l;
                    l = "-" == g.srcOri ? g.destMax - b : b - g.destMin;
                    for (var q = g.blocks, h = 0; h < q.length; ++h) {
                        var t = q[h],
                            c = t[0],
                            w = t[1],
                            t = t[2];
                        if (l >= w && l <= w + t) return d =
                            l - w, f = d + c + g.srcMin, f = "-" == g.destOri ? g.srcMax - c - d : d + c + g.srcMin, {
                                seq: g.srcChr,
                                pos: f,
                                flipped: g.srcOri != g.destOri
                            }
                    }
                }
            }
            return null
        };
        p.prototype.sourceBlocksForRange = function(a, d, f, g) {
            for (var l = this, h = d / this.granularity | 0, t = f / this.granularity | 0, e = !1, m = h; m <= t; ++m) {
                var c = a + "_" + m;
                2 != this.fetchedTiles[c] && (e = !0, 1 != this.fetchedTiles[c] && (this.fetchedTiles[c] = 1))
            }
            if (e) this.postFetchQueues[a] || this.chainFetcher.fetchChains(a, h * this.granularity, (t + 1) * this.granularity - 1).then(function(c) {
                l.chainsByDest || (l.chainsByDest[a] = []);
                for (var b = 0; b < c.length; ++b) {
                    var d = c[b],
                        k = l.chainsBySrc[d.srcChr];
                    if (k) {
                        for (var f = !1, g = 0; g < k.length; ++g) {
                            var w = k[g];
                            if (w.srcMin == d.srcMin && w.srcMax == d.srcMax) {
                                f = !0;
                                break
                            }
                        }
                        f || k.push(d)
                    } else l.chainsBySrc[d.srcChr] = [d];
                    if (k = l.chainsByDest[d.destChr]) {
                        f = !1;
                        for (g = 0; g < k.length; ++g)
                            if (w = k[g], w.destMin == d.destMin && w.destMax == d.destMax) {
                                f = !0;
                                break
                            }
                        f || k.push(d)
                    } else l.chainsByDest[d.destChr] = [d]
                }
                for (c = h; c <= t; ++c) l.fetchedTiles[a + "_" + c] = 2;
                if (l.postFetchQueues[a]) {
                    c = l.postFetchQueues[a];
                    for (b = 0; b < c.length; ++b) c[b]();
                    l.postFetchQueues[a] = null
                }
            }).catch(function(a) {
                console.log(a)
            }), q(this.postFetchQueues, a, function() {
                l.sourceBlocksForRange(a, d, f, g)
            });
            else {
                e = [];
                m = this.chainsByDest[a] || [];
                for (c = 0; c < m.length; ++c) {
                    var w = m[c];
                    if (d <= w.destMax && f >= w.destMin) {
                        var B, r;
                        "-" == w.srcOri ? (B = w.destMax - f, r = w.destMax - d) : (B = d - w.destMin, r = f - w.destMin);
                        for (var x = w.blocks, n = 0; n < x.length; ++n) {
                            var A = x[n],
                                p = A[0],
                                k = A[1],
                                E = A[2];
                            r >= k && B <= k + E && (A = Math.max(B, k) - k, k = Math.min(r, k + E) - k, "-" == w.destOri ? e.push(new b(w.srcChr, w.srcMax - p - k, w.srcMax -
                                p - A)) : e.push(new b(w.srcChr, w.srcMin + A + p, w.srcMin + k + p)))
                        }
                    }
                }
                g(e)
            }
        };
        n.prototype.fetchChains = function(a, b, d) {
            var f = this;
            return new t(function(b, d) {
                f.source.alignments(a, {}, function(a) {
                    for (var d = [], g = 0; g < a.length; ++g)
                        for (var c = a[g], w = 0; w < c.blocks.length; ++w) {
                            for (var q = c.blocks[w], h, t, e = 0; e < q.segments.length; ++e) {
                                var m = q.segments[e],
                                    r = c.objects[m.object];
                                r.dbSource === f.srcTag ? h = m : r.dbSource === f.destTag && (t = m)
                            }
                            if (h && t) {
                                for (var q = {
                                        srcChr: c.objects[h.object].accession,
                                        srcMin: h.min | 0,
                                        srcMax: h.max | 0,
                                        srcOri: h.strand,
                                        destChr: c.objects[t.object].accession,
                                        destMin: t.min | 0,
                                        destMax: t.max | 0,
                                        destOri: t.strand,
                                        blocks: []
                                    }, e = l(h.cigar), m = l(t.cigar), k = r = 0, E = 0, n = 0; E < e.length && n < m.length;)
                                    if ("M" == e[E].op && "M" == m[n].op) {
                                        var p = Math.min(e[E].cnt, m[n].cnt);
                                        q.blocks.push([r, k, p]);
                                        e[E].cnt == p ? ++E : e[E].cnt -= p;
                                        m[n].cnt == p ? ++n : m[n] -= p;
                                        r += p;
                                        k += p
                                    } else "I" == e[E].op ? k += e[E++].cnt : "I" == m[n].op && (r += m[n++].cnt);
                                d.push(q)
                            }
                        }
                    b(d)
                })
            })
        };
        m.prototype.fetchChains = function(a, b, d) {
            return this.bwg.then(function(f, g) {
                if (!f) throw Error("No BWG");
                return new t(function(g,
                    l) {
                    f.getUnzoomedView().readWigData(a, b, d, function(a) {
                        g(a.map(r))
                    })
                })
            })
        };
        z.prototype.fetchChains = function(a, b, d) {
            b = [];
            d = this.forwardAliases[a] || [];
            for (var f = 0; f < d.length; ++f) b.push({
                srcChr: d[f],
                srcMin: 1,
                srcMax: 1E9,
                srcOri: "+",
                destChr: a,
                destMin: 1,
                destMax: 1E9,
                destOri: "+",
                blocks: [
                    [1, 1, 1E9]
                ]
            });
            return t.resolve(b)
        };
        "undefined" !== typeof u && (u.exports = {
            Chainset: p
        })
    }, {
        "./bigwig": 3,
        "./bin": 4,
        "./cigar": 8,
        "./das": 10,
        "./utils": 48,
        "es6-promise": 52
    }],
    8: [function(e, u, s) {
        function p(e) {
            for (var h = [], p; null != (p = n.exec(e));) {
                var r =
                    p[1];
                0 == r.length && (r = 1);
                h.push({
                    cnt: r | 0,
                    op: p[2]
                })
            }
            return h
        }
        var n = /([0-9]*)([MIDS])/g;
        "undefined" !== typeof u && (u.exports = {
            parseCigar: p
        })
    }, {}],
    9: [function(e, u, s) {
        function p(a, q, g, l) {
            this.red = a | 0;
            this.green = q | 0;
            this.blue = g | 0;
            l && (this.name = l)
        }

        function n(a) {
            a = "00" + a.toString(16);
            return a.substring(a.length - 2)
        }

        function m(b) {
            var q = r[b];
            q || ((q = z.exec(b)) ? q = new p("0x" + q[1] | 0, "0x" + q[2] | 0, "0x" + q[3] | 0, b) : (q = a.exec(b)) ? q = new p(q[1] | 0, q[2] | 0, q[3] | 0, b) : (console.log("couldn't handle color: " + b), q = r.black), r[b] = q);
            return q
        }

        function h(a, q, g) {
            for (var l = [], f = 0; f < g.length; ++f) l.push(m(g[f]));
            g = [];
            f = 0;
            a: for (; f < a; ++f) {
                for (var d = q[0] + 1 * f / (a - 1) * (q[q.length - 1] - q[0]), h = 0; h < q.length - 1; ++h)
                    if (d >= q[h] && d <= q[h + 1]) {
                        var d = (d - q[h]) / (q[h + 1] - q[h]),
                            e = l[h],
                            h = l[h + 1],
                            h = (new p(e.red * (1 - d) + h.red * d | 0, e.green * (1 - d) + h.green * d | 0, e.blue * (1 - d) + h.blue * d | 0)).toSvgString();
                        g.push(h);
                        continue a
                    }
                throw "Bad step";
            }
            return g
        }

        function v(a, q, g, l) {
            return l ? h(a, [0, 0.5, 1], [q, g, l]) : h(a, [0, 1], [q, g])
        }
        p.prototype.toSvgString = function() {
            this.name || (this.name =
                "rgb(" + this.red + "," + this.green + "," + this.blue + ")");
            return this.name
        };
        p.prototype.toHexString = function() {
            return "#" + n(this.red) + n(this.green) + n(this.blue)
        };
        var r = {
                red: new p(255, 0, 0, "red"),
                green: new p(0, 255, 0, "green"),
                blue: new p(0, 0, 255, "blue"),
                yellow: new p(255, 255, 0, "yellow"),
                white: new p(255, 255, 255, "white"),
                black: new p(0, 0, 0, "black"),
                gray: new p(180, 180, 180, "gray"),
                grey: new p(180, 180, 180, "grey"),
                lightskyblue: new p(135, 206, 250, "lightskyblue"),
                lightsalmon: new p(255, 160, 122, "lightsalmon"),
                hotpink: new p(255,
                    105, 180, "hotpink")
            },
            z = /^#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$/,
            a = /rgb\(([0-9]+),([0-9]+),([0-9]+)\)/;
        "undefined" !== typeof u && (u.exports = {
            makeColourSteps: h,
            makeGradient: v,
            dasColourForName: m
        })
    }, {}],
    10: [function(e, u, s) {
        function p(a, b, c, d) {
            this.name = a;
            this.start = b;
            this.end = c;
            this.description = d
        }

        function n(a, b) {
            var c;
            "string" == typeof a ? (this.uri = a, c = b || {}) : c = a || {};
            for (var d in c) "function" != typeof c[d] && (this[d] = c[d]);
            this.coords || (this.coords = []);
            this.props || (this.props = {});
            (this.dasBaseURI =
                this.uri) && "/" != this.dasBaseURI.substr(this.uri.length - 1) && (this.dasBaseURI += "/")
        }

        function m() {}

        function h(a, b) {
            return a.taxon == b.taxon && a.auth == b.auth && a.version == b.version
        }

        function v(a, b, c, d, f) {
            this.name = a;
            this.start = b;
            this.end = c;
            this.alphabet = d;
            this.seq = f
        }

        function r() {}

        function z(a) {
            a && (this.id = a)
        }

        function a(a, b) {
            this.desc = a;
            this.uri = b
        }

        function b(a) {
            this.type = a;
            this.objects = {};
            this.blocks = []
        }

        function q() {
            this.styles = []
        }

        function g() {}

        function l(a, b) {
            b = b || {};
            this.uri = a;
            this.opts = b
        }

        function f(a, b) {
            var c =
                a.getElementsByTagName(b);
            if (0 < c.length && c[0].firstChild) {
                c = c[0];
                if (1 == c.childNodes.length) return c.firstChild.nodeValue;
                for (var d = "", f = 0; f < c.childNodes.length; ++f) d += c.childNodes[f].nodeValue;
                return d
            }
            return null
        }

        function d(b) {
            for (var d = [], c = b.getElementsByTagName("LINK"), f = 0; f < c.length; ++f) {
                var g = c[f];
                g.parentNode == b && d.push(new a(g.firstChild ? g.firstChild.nodeValue : "Unknown", g.getAttribute("href")))
            }
            return d
        }

        function t(a) {
            var b = [];
            a = a.getElementsByTagName("NOTE");
            for (var c = 0; c < a.length; ++c) a[c].firstChild &&
                b.push(a[c].firstChild.nodeValue);
            return b
        }

        function D(a, b, c, d) {
            if (window.XDomainRequest) {
                var f = new XDomainRequest;
                f.onload = function() {
                    var a = new ActiveXObject("Microsoft.XMLDOM");
                    a.async = !1;
                    a.loadXML(f.responseText);
                    b(a)
                };
                f.open("get", a)
            } else Date.now(), f = new XMLHttpRequest, f.onreadystatechange = function() {
                4 == f.readyState && (200 <= f.status || 0 == f.status) && b(f.responseXML, f)
            }, f.open("get", a, !0), c && (f.withCredentials = !0), d && f.setRequestHeader("X-DAS-Authorisation", d), f.overrideMimeType("text/xml"), f.setRequestHeader("Accept",
                "application/xml,*/*");
            f.send("")
        }

        function P(a) {
            a = ("" + a).toLowerCase();
            return "yes" === a || "true" === a
        }

        function G(a) {
            if (!a) return !1;
            a = ("" + a).toLowerCase();
            return "no" !== a || "false" !== a
        }

        function L(a) {
            var b = N(a);
            b.styles = [];
            for (var c = 0; c < a.styles.length; ++c) {
                var d = b.styles[c] = N(a.styles[c]);
                d._methodRE = d._labelRE = d._typeRE = void 0;
                d.style = N(d.style);
                d.style.id = void 0;
                d.style._gradient = void 0
            }
            return b
        }
        if ("undefined" !== typeof e) {
            s = e("./utils");
            var N = s.shallowCopy,
                Q = s.pusho,
                y = e("./color").makeColourSteps
        }
        p.prototype.toString =
            function() {
                return this.name + ":" + this.start + ".." + this.end
            };
        p.prototype.isBounded = function() {
            return this.start && this.end
        };
        p.prototype.toDASQuery = function() {
            var a = "segment=" + this.name;
            this.start && this.end && (a += ":" + this.start + "," + this.end);
            return a
        };
        n.prototype.entryPoints = function(a) {
            this.doCrossDomainRequest(this.dasBaseURI + "entry_points", function(b) {
                if (!b) return a([]);
                var c = [];
                b = b.getElementsByTagName("SEGMENT");
                for (var d = 0; d < b.length; ++d) {
                    var f = b[d],
                        g = f.getAttribute("id"),
                        l = f.getAttribute("size"),
                        q;
                    l ? (q = 1, l |= 0) : ((q = f.getAttribute("start")) && (q |= 0), (l = f.getAttribute("stop")) && (l |= 0));
                    var h = null;
                    f.firstChild && (h = f.firstChild.nodeValue);
                    c.push(new p(g, q, l, h))
                }
                a(c)
            })
        };
        n.prototype.sequence = function(a, b) {
            var c = this.dasBaseURI + "sequence?" + a.toDASQuery();
            this.doCrossDomainRequest(c, function(a) {
                if (a) {
                    var c = [];
                    a = a.getElementsByTagName("SEQUENCE");
                    for (var d = 0; d < a.length; ++d) {
                        var f = a[d],
                            g = f.getAttribute("id"),
                            l = f.getAttribute("start"),
                            q = f.getAttribute("stop"),
                            k = null;
                        if (f.firstChild)
                            for (var f = f.firstChild.nodeValue,
                                    k = "", h = 0;;) {
                                var t = f.indexOf("\n", h);
                                if (0 <= t) k += f.substring(h, t).toUpperCase(), h = t + 1;
                                else {
                                    k += f.substring(h).toUpperCase();
                                    break
                                }
                            }
                        c.push(new v(g, l, q, "DNA", k))
                    }
                    b(c)
                } else b([])
            })
        };
        n.prototype.features = function(a, b, c) {
            b = b || {};
            var g;
            if (this.features_uri) g = this.features_uri;
            else {
                var l = [];
                if (a) l.push(a.toDASQuery());
                else if (b.group)
                    if (a = b.group, "string" == typeof a) l.push("group_id=" + a);
                    else
                        for (var q = 0; q < a.length; ++q) l.push("group_id=" + a[q]);
                if (b.adjacent)
                    for (a = b.adjacent, "string" == typeof a && (a = [a]), q = 0; q <
                        a.length; ++q) l.push("adjacent=" + a[q]);
                if (b.type)
                    if ("string" == typeof b.type) l.push("type=" + b.type);
                    else
                        for (a = 0; a < b.type.length; ++a) l.push("type=" + b.type[a]);
                b.maxbins && l.push("maxbins=" + b.maxbins);
                0 < l.length ? g = this.dasBaseURI + "features?" + l.join(";") : c([], "No filters specified")
            }
            this.doCrossDomainRequest(g, function(a, b) {
                if (a) {
                    for (var g = [], l = {}, k = a.getElementsByTagName("SEGMENT"), q = 0; q < k.length; ++q) {
                        var h = k[q],
                            w = h.getAttribute("id");
                        l[w] = {
                            min: h.getAttribute("start"),
                            max: h.getAttribute("stop")
                        };
                        for (var h =
                                h.getElementsByTagName("FEATURE"), e = 0; e < h.length; ++e) {
                            var m = h[e],
                                B = new r;
                            B.segment = w;
                            B.id = m.getAttribute("id");
                            B.label = m.getAttribute("label");
                            var n = f(m, "START"),
                                p = f(m, "END");
                            (n | 0) > (p | 0) ? (B.min = p | 0, B.max = n | 0) : (B.min = n | 0, B.max = p | 0);
                            n = m.getElementsByTagName("TYPE");
                            0 < n.length && (n = n[0], n.firstChild && (B.type = n.firstChild.nodeValue), B.typeId = n.getAttribute("id"), B.typeCv = n.getAttribute("cvId"));
                            B.type = f(m, "TYPE");
                            !B.type && B.typeId && (B.type = B.typeId);
                            B.method = f(m, "METHOD");
                            (n = f(m, "ORIENTATION")) || (n = "0");
                            B.orientation = n;
                            B.score = f(m, "SCORE");
                            B.links = d(m);
                            B.notes = t(m);
                            n = m.getElementsByTagName("GROUP");
                            for (p = 0; p < n.length; ++p) {
                                var y = n[p],
                                    D = new z;
                                D.type = y.getAttribute("type");
                                D.id = y.getAttribute("id");
                                D.links = d(y);
                                D.notes = t(y);
                                B.groups ? B.groups.push(D) : B.groups = Array(D)
                            }
                            if (B.notes)
                                for (n = 0; n < B.notes.length; ++n) p = B.notes[n], 0 == p.indexOf("Genename=") && (y = new z, y.type = "gene", y.id = p.substring(9), B.groups ? B.groups.push(y) : B.groups = Array(y));
                            n = m.getElementsByTagName("PART");
                            if (0 < n.length) {
                                y = [];
                                for (p = 0; p < n.length; ++p) y.push(n[p].getAttribute("id"));
                                B.parts = y
                            }
                            n = m.getElementsByTagName("PARENT");
                            if (0 < n.length) {
                                m = [];
                                for (p = 0; p < n.length; ++p) m.push(n[p].getAttribute("id"));
                                B.parents = m
                            }
                            g.push(B)
                        }
                    }
                    c(g, void 0, l)
                } else c([], "Failed request: " + (0 == b.status ? "server may not support CORS" : "status=" + b.status))
            }, function(a) {
                c([], a)
            })
        };
        n.prototype.alignments = function(a, d, c) {
            var g = this.dasBaseURI + "alignment?query=" + a;
            this.doCrossDomainRequest(g, function(a) {
                if (a) {
                    var d = [];
                    a = a.getElementsByTagName("alignment");
                    for (var l = 0; l < a.length; ++l) {
                        for (var q = a[l], h = new b(q.getAttribute("alignType")),
                                t = q.getElementsByTagName("alignObject"), k = 0; k < t.length; ++k) {
                            var e = t[k],
                                e = {
                                    id: e.getAttribute("intObjectId"),
                                    accession: e.getAttribute("dbAccessionId"),
                                    version: e.getAttribute("objectVersion"),
                                    dbSource: e.getAttribute("dbSource"),
                                    dbVersion: e.getAttribute("dbVersion")
                                };
                            h.objects[e.id] = e
                        }
                        q = q.getElementsByTagName("block");
                        for (t = 0; t < q.length; ++t) {
                            for (var e = q[t], k = {
                                    order: e.getAttribute("blockOrder"),
                                    segments: []
                                }, e = e.getElementsByTagName("segment"), m = 0; m < e.length; ++m) {
                                var r = e[m],
                                    r = {
                                        object: r.getAttribute("intObjectId"),
                                        min: r.getAttribute("start"),
                                        max: r.getAttribute("end"),
                                        strand: r.getAttribute("strand"),
                                        cigar: f(r, "cigar")
                                    };
                                k.segments.push(r)
                            }
                            h.blocks.push(k)
                        }
                        d.push(h)
                    }
                    c(d)
                } else c([], "Failed request " + g)
            })
        };
        q.prototype.pushStyle = function(a, b, c) {
            a || (a = {
                type: "default"
            });
            a = N(a);
            b && (a.zoom = b);
            a.style = c;
            this.styles.push(a)
        };
        n.prototype.stylesheet = function(a, b) {
            var c, d = this.credentials;
            this.stylesheet_uri ? (c = this.stylesheet_uri, d = !1) : c = this.dasBaseURI + "stylesheet";
            D(c, function(c) {
                if (c) {
                    var d = new q;
                    c = c.getElementsByTagName("TYPE");
                    for (var f = 0; f < c.length; ++f) {
                        var l = c[f],
                            h = {};
                        h.type = l.getAttribute("id");
                        h.label = l.getAttribute("label");
                        h.method = l.getAttribute("method");
                        for (var l = l.getElementsByTagName("GLYPH"), w = 0; w < l.length; ++w) {
                            var k = l[w],
                                t = k.getAttribute("zoom"),
                                e;
                            a: {
                                if (k.hasChildNodes()) {
                                    k = k.firstChild;
                                    do {
                                        if (k.nodeType == Node.ELEMENT_NODE) {
                                            e = k;
                                            break a
                                        }
                                        k = k.nextSibling
                                    } while (null != k)
                                }
                                e = null
                            }
                            k = new g;
                            k.glyph = e.localName;
                            for (e = e.firstChild; e;) {
                                if (e.nodeType == Node.ELEMENT_NODE)
                                    if ("BGGRAD" == e.localName) {
                                        var m = k,
                                            r = e.localName,
                                            n, p =
                                            e;
                                        n = (n = p.getAttribute("steps")) ? n | 0 : 50;
                                        for (var D = [], z = [], p = p.getElementsByTagName("STOP"), L = 0; L < p.length; ++L) {
                                            var v = p[L];
                                            D.push(1 * v.getAttribute("score"));
                                            z.push(v.firstChild.nodeValue)
                                        }
                                        n = y(n, D, z);
                                        m[r] = n
                                    } else k[e.localName] = e.firstChild.nodeValue;
                                e = e.nextSibling
                            }
                            d.pushStyle(h, t, k)
                        }
                    }
                    a(d)
                } else b && b()
            }, d)
        };
        l.prototype.sources = function(a, b, c) {
            c || (c = {});
            var d = [];
            c.taxon && d.push("organism=" + c.taxon);
            c.auth && d.push("authority=" + c.auth);
            c.version && d.push("version=" + c.version);
            c = this.uri;
            0 < d.length && (c = c +
                "?" + d.join("&"));
            D(c, function(c) {
                if (!c && b) b();
                else {
                    var d = [];
                    c = c.getElementsByTagName("SOURCE");
                    for (var f = 0; f < c.length; ++f) {
                        var g = c[f],
                            l = g.getElementsByTagName("VERSION");
                        if (!(1 > l.length)) {
                            for (var q = l[0], k = q.getElementsByTagName("COORDINATES"), l = [], h = 0; h < k.length; ++h) {
                                var w = k[h],
                                    e = new m;
                                e.auth = w.getAttribute("authority");
                                e.taxon = w.getAttribute("taxid");
                                e.version = w.getAttribute("version");
                                l.push(e)
                            }
                            for (var k = [], w = q.getElementsByTagName("CAPABILITY"), t, h = 0; h < w.length; ++h) e = w[h], k.push(e.getAttribute("type")),
                                "das1:features" == e.getAttribute("type") && (t = e.getAttribute("query_uri"), t = t.substring(0, t.length - 8));
                            h = {};
                            q = q.getElementsByTagName("PROP");
                            for (w = 0; w < q.length; ++w) Q(h, q[w].getAttribute("name"), q[w].getAttribute("value"));
                            t && (g = new n(t, {
                                source_uri: g.getAttribute("uri"),
                                name: g.getAttribute("title"),
                                desc: g.getAttribute("description"),
                                coords: l,
                                props: h,
                                capabilities: k
                            }), d.push(g))
                        }
                    }
                    a(d)
                }
            })
        };
        n.prototype.doCrossDomainRequest = function(a, b, c) {
            var d;
            this.xUser && (d = "Basic " + btoa(this.xUser + ":" + this.xPass));
            try {
                return D(a,
                    b, this.credentials, d)
            } catch (f) {
                if (c) c(f);
                else throw f;
            }
        };
        "undefined" !== typeof u && (u.exports = {
            DASGroup: z,
            DASFeature: r,
            DASStylesheet: q,
            DASStyle: g,
            DASSource: n,
            DASSegment: p,
            DASRegistry: l,
            DASSequence: v,
            DASLink: a,
            isDasBooleanTrue: P,
            isDasBooleanNotFalse: G,
            copyStylesheet: L,
            coordsMatch: h
        })
    }, {
        "./color": 9,
        "./utils": 48
    }],
    11: [function(e, u, s) {
        function p(h, e, n) {
            function a() {
                n ? (b.className = "fa fa-caret-down", e.style.display = "table") : (b.className = "fa fa-caret-right", e.style.display = "none")
            }
            var b = m("i");
            a();
            b.addEventListener("click",
                function(b) {
                    b.preventDefault();
                    b.stopPropagation();
                    n = !n;
                    a()
                }, !1);
            h = m("h6", [b, " ", h], {}, {
                display: "block",
                background: "gray",
                color: "white",
                width: "100%",
                padding: "5px 2px",
                margin: "0px"
            });
            return m("div", [h, e], {})
        }
        if ("undefined" !== typeof e) {
            var n = e("./cbrowser").Browser;
            e = e("./utils");
            var m = e.makeElement,
                h = e.removeChildren
        }
        n.prototype.removeAllPopups = function() {
            h(this.hPopupHolder);
            h(this.popupHolder)
        };
        n.prototype.makeTooltip = function(h, e) {
            var n = !1,
                a = this,
                b = null,
                q;
            q = function(a) {
                n = !1;
                b && (clearTimeout(b),
                    b = null);
                h.removeEventListener("mouseout", q, !1)
            };
            var g = function(l) {
                var f = l.clientX + window.scrollX,
                    d = l.clientY + window.scrollY;
                b || (b = setTimeout(function() {
                    var l;
                    l = "function" === typeof e ? e() : e;
                    var q = m("div", [m("div", null, {
                        className: "tooltip-arrow"
                    }), m("div", l, {
                        className: "tooltip-inner"
                    })], {
                        className: "tooltip bottom in"
                    }, {
                        display: "block",
                        top: "" + (d + 20) + "px",
                        left: "" + Math.max(f - 30, 20) + "px"
                    });
                    a.hPopupHolder.appendChild(q);
                    var p;
                    p = function(b) {
                        try {
                            a.hPopupHolder.removeChild(q)
                        } catch (d) {}
                        window.removeEventListener("mousemove",
                            p, !1);
                        n && null != h.offsetParent && g(b)
                    };
                    window.addEventListener("mousemove", p, !1);
                    b = null
                }, 1E3))
            };
            h.addEventListener("mouseover", function(a) {
                n = !0;
                h.addEventListener("mouseout", q, !1);
                g(a)
            }, !1);
            h.addEventListener("DOMNodeRemovedFromDocument", function(a) {
                n = !1;
                b && (clearTimeout(b), b = null)
            }, !1)
        };
        n.prototype.popit = function(h, e, n, a) {
            var b = this;
            a || (a = {});
            h || (h = {});
            var q = a.width || 200,
                g;
            h.clientX ? (a = h.clientX, g = h.clientY) : (a = 500, g = 50);
            a += document.documentElement.scrollLeft || document.body.scrollLeft;
            g += document.documentElement.scrollTop ||
                document.body.scrollTop;
            var l = window.innerWidth,
                f = g,
                d = Math.min(a - q / 2 - 4, l - q - 30),
                t = m("div");
            t.className = "popover fade " + (h.clientX ? "bottom " : "") + "in";
            t.style.display = "block";
            t.style.position = "absolute";
            t.style.top = "" + f + "px";
            t.style.left = "" + d + "px";
            t.style.width = q + "px";
            276 < q && (t.style.maxWidth = q + "px");
            t.appendChild(m("div", null, {
                className: "arrow"
            }));
            if (e) {
                var p = m("button", "", {
                    className: "close"
                });
                p.innerHTML = "&times;";
                p.addEventListener("mouseover", function(a) {
                    p.style.color = "red"
                }, !1);
                p.addEventListener("mouseout",
                    function(a) {
                        p.style.color = "black"
                    }, !1);
                p.addEventListener("click", function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    b.removeAllPopups()
                }, !1);
                h = m("h4", [m("span", e, null, {
                    maxWidth: "200px"
                }), p], {}, {
                    paddingLeft: "10px",
                    paddingRight: "10px"
                });
                var P, G, L, N;
                L = function(a) {
                    a.stopPropagation();
                    a.preventDefault();
                    d += a.clientX - P;
                    8 > d && (d = 8);
                    d > l - q - 32 && (d = l - q - 26);
                    f += a.clientY - G;
                    f = Math.max(10, f);
                    t.style.top = "" + f + "px";
                    t.style.left = "" + Math.min(d, l - q - 10) + "px";
                    P = a.clientX;
                    G = a.clientY
                };
                N = function(a) {
                    a.stopPropagation();
                    a.preventDefault();
                    window.removeEventListener("mousemove", L, !1);
                    window.removeEventListener("mouseup", N, !1)
                };
                h.addEventListener("mousedown", function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    P = a.clientX;
                    G = a.clientY;
                    window.addEventListener("mousemove", L, !1);
                    window.addEventListener("mouseup", N, !1)
                }, !1);
                t.appendChild(h)
            }
            t.appendChild(m("div", n, {
                className: "popover-content"
            }, {
                padding: "0px"
            }));
            this.hPopupHolder.appendChild(t);
            var Q = {
                node: t,
                displayed: !0
            };
            t.addEventListener("DOMNodeRemoved", function(a) {
                a.target ==
                    t && (Q.displayed = !1)
            }, !1);
            return Q
        };
        "undefined" !== typeof u && (u.exports = {
            makeTreeTableSection: p
        })
    }, {
        "./cbrowser": 6,
        "./utils": 48
    }],
    12: [function(e, u, s) {
        function p(e, m) {
            0 > e.indexOf("?") && (e += "?soft=true");
            return new h(function(h, a) {
                var b = new XMLHttpRequest;
                b.onreadystatechange = function() {
                    if (4 == b.readyState)
                        if (300 <= b.status) a("Error code " + b.status);
                        else {
                            var q = JSON.parse(b.response);
                            h(m ? q : q.location)
                        }
                };
                b.open("GET", e, !0);
                b.setRequestHeader("Accept", "application/json");
                b.responseType = "text";
                b.send("")
            })
        }

        function n(h) {
            this.rawurl = h
        }

        function m(h, e, m, a) {
            a || ("object" === typeof e ? (a = e, e = void 0) : a = {});
            this.url = "string" === typeof h ? new n(h) : h;
            this.start = e || 0;
            m && (this.end = m);
            this.opts = a
        }
        if ("undefined" !== typeof e) var h = e("es6-promise").Promise;
        n.prototype.getURLPromise = function() {
            this.urlPromise && this.urlPromiseValidity > Date.now() || (this.urlPromise = p(this.rawurl, !0).then(function(h) {
                return h.location
            }), this.urlPromiseValidity = Date.now() + 432E5);
            return this.urlPromise
        };
        m.prototype.slice = function(h, e) {
            if (0 > h) throw "Bad slice " +
                h;
            var n = this.start,
                a = this.end,
                n = n && h ? n + h : h || n;
            return new m(this.url, n, e && n ? n + e - 1 : a || e - 1, this.opts)
        };
        m.prototype.fetchAsText = function(h) {
            var e = this,
                m = new XMLHttpRequest;
            e.url.getURLPromise().then(function(a) {
                m.open("GET", a, !0);
                if (e.end) {
                    if (1E8 < e.end - e.start) throw "Monster fetch!";
                    m.setRequestHeader("Range", "bytes=" + e.start + "-" + e.end)
                }
                m.onreadystatechange = function() {
                    if (4 == m.readyState) return 200 == m.status || 206 == m.status ? h(m.responseText) : h(null)
                };
                e.opts.credentials && (m.withCredentials = !0);
                m.send("")
            }).catch(function(a) {
                console.log(a);
                return h(null)
            })
        };
        m.prototype.salted = function() {
            return this
        };
        m.prototype.fetch = function(h, e, m) {
            var a = this;
            e = e || 1;
            if (3 < e) return h(null);
            a.url.getURLPromise().then(function(b) {
                var q = new XMLHttpRequest,
                    g;
                q.open("GET", b, !0);
                q.overrideMimeType("text/plain; charset=x-user-defined");
                if (a.end) {
                    if (1E8 < a.end - a.start) throw "Monster fetch!";
                    q.setRequestHeader("Range", "bytes=" + a.start + "-" + a.end);
                    g = a.end - a.start + 1
                }
                q.responseType = "arraybuffer";
                q.onreadystatechange = function() {
                    if (4 == q.readyState) {
                        if (200 == q.status ||
                            206 == q.status) {
                            if (q.response) {
                                var b = q.response.byteLength;
                                return !g || g == b || m && b == m ? h(q.response) : a.fetch(h, e + 1, b)
                            }
                            if (q.mozResponseArrayBuffer) return h(q.mozResponseArrayBuffer);
                            b = q.responseText;
                            if (!g || g == b.length || m && b.length == m) {
                                if (b = q.responseText) {
                                    for (var f = new Uint8Array(b.length), d = 0; d < f.length; ++d) f[d] = b.charCodeAt(d);
                                    b = f.buffer
                                } else b = null;
                                return h(b)
                            }
                            return a.fetch(h, e + 1, b.length)
                        }
                        return a.fetch(h, e + 1)
                    }
                };
                a.opts.credentials && (q.withCredentials = !0);
                q.send("")
            }).catch(function(a) {
                console.log(a)
            })
        };
        "undefined" !== typeof u && (u.exports = {
            lookupEncodeURI: p,
            EncodeFetchable: m
        })
    }, {
        "es6-promise": 52
    }],
    13: [function(e, u, s) {
        function p(h) {
            this.source = h;
            this.base = h.uri || "//rest.ensembl.org";
            if (0 === this.base.indexOf("//")) {
                var a = window.location.protocol;
                "http:" != a && "https:" != a && (this.base = "http:" + this.base)
            }
            this.species = h.species || "human";
            this.activityListeners = [];
            this.busy = 0;
            this.type = "string" === typeof h.type ? [h.type] : h.type || ["regulatory"]
        }
        if ("undefined" !== typeof e) {
            var n = e("./sourceadapters").registerSourceAdapterFactory;
            e = e("./das");
            var m = e.DASStylesheet,
                h = e.DASStyle,
                v = e.DASFeature,
                r = e.DASGroup
        }
        p.prototype.addActivityListener = function(h) {
            this.activityListeners.push(h)
        };
        p.prototype.notifyActivity = function() {
            for (var h = 0; h < this.activityListeners.length; ++h) try {
                this.activityListeners[h](this.busy)
            } catch (a) {
                console.log(a)
            }
        };
        p.prototype.getStyleSheet = function(e) {
            var a = new m,
                b = new h;
            b.glyph = "__NONE";
            0 <= this.type.indexOf("exon") && a.pushStyle({
                type: "transcript"
            }, null, b);
            (0 <= this.type.indexOf("exon") || 0 <= this.type.indexOf("transcript")) &&
            a.pushStyle({
                type: "gene"
            }, null, b);
            b = new h;
            b.glyph = "BOX";
            b.FGCOLOR = "black";
            b.BGCOLOR = "red";
            b.HEIGHT = 8;
            b.BUMP = !0;
            b.LABEL = !0;
            b.ZINDEX = 10;
            a.pushStyle({
                type: "cds"
            }, null, b);
            b = new h;
            b.glyph = "SQUARE";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.FGCOLOR = "blue";
            a.pushStyle({
                type: "variation",
                method: ".+_UTR_variant"
            }, null, b);
            b = new h;
            b.glyph = "TRIANGLE";
            b.DIRECTION = "S";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.FGCOLOR = "blue";
            a.pushStyle({
                type: "variation",
                method: "missense_variant"
            }, null, b);
            b = new h;
            b.glyph = "TRIANGLE";
            b.DIRECTION = "N";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.FGCOLOR = "blue";
            a.pushStyle({
                type: "variation",
                method: "splice_.+_variant"
            }, null, b);
            b = new h;
            b.glyph = "STAR";
            b.POINTS = 6;
            b.BUMP = "yes";
            b.LABEL = "no";
            b.FGCOLOR = "blue";
            a.pushStyle({
                type: "variation",
                method: "regulatory_region_variant"
            }, null, b);
            b = new h;
            b.glyph = "PLIMSOLL";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.FGCOLOR = "rgb(50,80,255)";
            b.STROKECOLOR = "black";
            a.pushStyle({
                type: "variation"
            }, null, b);
            b = new h;
            b.glyph = "SQUARE";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.BGCOLOR = "#888888";
            b.FGCOLOR = "red";
            a.pushStyle({
                type: "indel",
                method: ".+_UTR_variant"
            }, null, b);
            b = new h;
            b.glyph = "TRIANGLE";
            b.DIRECTION = "S";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.BGCOLOR = "#888888";
            b.FGCOLOR = "red";
            a.pushStyle({
                type: "indel",
                method: "missense_variant"
            }, null, b);
            b = new h;
            b.glyph = "TRIANGLE";
            b.DIRECTION = "N";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.BGCOLOR = "#888888";
            b.FGCOLOR = "red";
            a.pushStyle({
                type: "indel",
                method: "splice_.+_variant"
            }, null, b);
            b = new h;
            b.glyph = "STAR";
            b.POINTS = 6;
            b.BUMP = "yes";
            b.LABEL = "no";
            b.BGCOLOR = "#888888";
            b.FGCOLOR = "red";
            a.pushStyle({
                    type: "indel",
                    method: "regulatory_region_variant"
                },
                null, b);
            b = new h;
            b.glyph = "PLIMSOLL";
            b.BUMP = "yes";
            b.LABEL = "no";
            b.BGCOLOR = "#888888";
            b.FGCOLOR = "red";
            b.STROKECOLOR = "black";
            a.pushStyle({
                type: "indel"
            }, null, b);
            b = new h;
            b.glyph = "BOX";
            b.FGCOLOR = "black";
            b.BGCOLOR = "orange";
            b.HEIGHT = 8;
            b.BUMP = !0;
            b.LABEL = !0;
            b.ZINDEX = 20;
            a.pushStyle({
                type: "default"
            }, null, b);
            return e(a)
        };
        p.prototype.getScales = function() {
            return []
        };
        p.prototype.fetch = function(h, a, b, q, g, l, f) {
            var d = this;
            a = this.base + "/overlap/region/" + this.species + "/" + h + ":" + a + "-" + b;
            b = [];
            for (q = 0; q < this.type.length; ++q) b.push("feature=" +
                this.type[q]);
            b.push("content-type=application/json");
            a = a + "?" + b.join(";");
            var e = new XMLHttpRequest;
            e.onreadystatechange = function() {
                if (4 == e.readyState)
                    if (d.busy--, d.notifyActivity(), 300 <= e.status) {
                        var a = "Error code " + e.status;
                        try {
                            var b = JSON.parse(e.response);
                            b.error && (a = b.error)
                        } catch (g) {}
                        f(a, null)
                    } else {
                        for (var a = JSON.parse(e.response), b = [], l = 0; l < a.length; ++l) {
                            var q = a[l],
                                m = [],
                                n = new v;
                            n.segment = h;
                            n.min = q.start | 0;
                            n.max = q.end | 0;
                            n.type = q.feature_type || "unknown";
                            n.id = q.ID;
                            if (q.Parent) {
                                var p = new r;
                                p.id =
                                    q.Parent;
                                n.groups = [p]
                            }
                            q.strand && (0 > q.strand ? n.orientation = "-" : 0 < q.strand && (n.orientation = "+"));
                            q.consequence_type && (n.method = q.consequence_type);
                            q.alt_alleles && (m.push("Alleles=" + q.alt_alleles.join("/")), 1 < q.alt_alleles.length && (q.alt_alleles[1].length != q.alt_alleles[0].length || "-" == q.alt_alleles[1]) && (n.type = "indel"));
                            0 < m.length && (n.notes = m);
                            b.push(n)
                        }
                        f(null, b)
                    }
            };
            d.busy++;
            d.notifyActivity();
            e.open("GET", a, !0);
            e.responseType = "text";
            e.send("")
        };
        n("ensembl", function(h) {
            return {
                features: new p(h)
            }
        })
    }, {
        "./das": 10,
        "./sourceadapters": 34
    }],
    14: [function(e, u, s) {
        if ("undefined" !== typeof e) var p = e("./cbrowser").Browser,
            n = e("./utils").shallowCopy,
            m = e("./sha1").hex_sha1,
            h = e("./das").copyStylesheet;
        p.prototype.exportFullConfig = function(h) {
            h = {
                chr: this.chr,
                viewStart: this.viewStart | 0,
                viewEnd: this.viewEnd | 0,
                cookieKey: "dalliance_" + m(Date.now()),
                coordSystem: this.coordSystem,
                sources: this.exportSourceConfig(),
                chains: this.exportChains()
            };
            this.prefix && (h.prefix = this.prefix);
            return h
        };
        p.prototype.exportChains = function() {
            var h = {},
                e = this.chains || {},
                m;
            for (m in e) h[m] = e[m].exportConfig();
            return h
        };
        p.prototype.exportSourceConfig = function(e) {
            e = [];
            for (var m = 0; m < this.tiers.length; ++m) {
                var p = this.tiers[m],
                    a = n(p.dasSource);
                a.noPersist || (a.coords = void 0, a.props = void 0, a.disabled || (a.disabled = void 0), p.config.stylesheet ? (a.style = h(p.config.stylesheet).styles, a.stylesheet_uri = void 0) : a.style && (a.style = h({
                        styles: a.style
                    }).styles), "string" === typeof p.config.name && (a.name = p.config.name), void 0 !== p.config.height && (a.forceHeight = p.config.height),
                    void 0 !== p.config.forceMin && (a.forceMin = p.config.forceMin), p.config.forceMinDynamic && (a.forceMinDynamic = p.config.forceMinDynamic), void 0 !== p.config.forceMax && (a.forceMax = p.config.forceMax), void 0 !== p.config.bumped && (a.bumped = p.config.bumped), p.config.forceMaxDynamic && (a.forceMaxDynamic = p.config.forceMaxDynamic), e.push(a))
            }
            return e
        };
        p.prototype.exportPageTemplate = function(h) {
            h = h || {};
            return '<html>\n  <head>\n    <script language="javascript" src="' + this.resolveURL("$$dalliance-compiled.js") + '">\x3c/script>\n    <script language="javascript">\n      var dalliance_browser = new Browser(' +
                JSON.stringify(this.exportFullConfig(h), null, 2) + ');\n    \x3c/script>\n  </head>\n  <body>\n    <div id="svgHolder">Dalliance goes here</div>\n  </body>\n<html>\n'
        }
    }, {
        "./cbrowser": 6,
        "./das": 10,
        "./sha1": 33,
        "./utils": 48
    }],
    15: [function(e, u, s) {
        if ("undefined" !== typeof e) var p = e("./cbrowser").Browser,
            n = e("./glyphs").OverlayLabelCanvas,
            m = e("./numformats").formatQuantLabel,
            h = e("./sequence-draw").drawSeqTierGC;
        p.prototype.exportImage = function(e) {
            e = e || {};
            for (var p = this.featurePanelWidth, z = 0, a = 0; a < this.tiers.length; ++a) {
                0 <
                    a && (z += 3);
                var b = this.tiers[a];
                void 0 !== b.layoutHeight && (z += b.layoutHeight)
            }
            var q = e.resolutionMultiplier || 1,
                a = (p + 200) * q | 0,
                b = z * q | 0,
                z = makeElement("canvas", null, {
                    width: a,
                    height: b
                }),
                g = z.getContext("2d");
            g.fillStyle = "white";
            g.fillRect(0, 0, a, b);
            g.scale(q, q);
            for (a = q = 0; a < this.tiers.length; ++a) {
                var b = this.tiers[a],
                    l = (b.glyphCacheOrigin - this.viewStart) * this.scale,
                    f = new n;
                g.save();
                g.translate(0, q);
                g.save();
                g.beginPath();
                g.moveTo(200, 0);
                g.lineTo(200 + p, 0);
                g.lineTo(200 + p, b.layoutHeight);
                g.lineTo(200, b.layoutHeight);
                g.closePath();
                g.clip();
                g.translate(200, 0);
                g.save();
                g.translate(l, 0);
                b.subtiers ? b.paintToContext(g, f, l + 1E3) : h(b, b.currentSequence, g);
                g.restore();
                g.save();
                g.translate(l, 0);
                f.draw(g, -l, p - l);
                g.restore();
                g.restore();
                for (var f = !1, d = 0, t = b.subtiers || [], D = 0; D < t.length; ++D) {
                    var P = t[D];
                    if (P.quant) {
                        var f = !0,
                            G = P.quant,
                            L = P.height,
                            N = 2;
                        40 < L && (N = 1 + (L / 20 | 0));
                        var L = L / (N - 1),
                            Q = (G.max - G.min) / (N - 1);
                        g.beginPath();
                        g.moveTo(205, d);
                        g.lineTo(200, d);
                        g.lineTo(200, d + P.height);
                        g.lineTo(205, d + P.height);
                        for (var y = 1; y < N - 1; ++y) {
                            var H =
                                y * L;
                            g.moveTo(200, d + H);
                            g.lineTo(203, d + H)
                        }
                        g.strokeStyle = "black";
                        g.strokeWidth = 2;
                        g.stroke();
                        g.fillStyle = "black";
                        var y = m(G.max),
                            F = d + 7;
                        g.fillText(y, 197 - g.measureText(y).width, F);
                        y = m(G.min);
                        F = d + P.height;
                        g.fillText(y, 197 - g.measureText(y).width, F);
                        for (y = 1; y < N - 1; ++y) H = y * L, F = m(1 * G.max - y * Q), H = d + H + 3, g.fillText(F, 197 - g.measureText(F).width, H)
                    }
                    d += P.height + 3
                }
                d = "string" === typeof b.config.name ? b.config.name : b.dasSource.name;
                t = g.measureText(d).width;
                g.fillStyle = "black";
                g.fillText(d, 200 - (f ? 22 : 12) - t, (b.layoutHeight +
                    6) / 2);
                g.restore();
                q += b.layoutHeight + 3
            }
            if (e.highlights) {
                g.save();
                g.beginPath();
                g.moveTo(200, 0);
                g.lineTo(200 + p, 0);
                g.lineTo(200 + p, q);
                g.lineTo(200, q);
                g.closePath();
                g.clip();
                g.translate(200 + l, 0);
                l = p = this.viewStart;
                a = this.viewEnd;
                for (b = 0; b < this.highlights.length; ++b) L = this.highlights[b], (L.chr === this.chr || L.chr === "chr" + this.chr) && L.min < a && L.max > l && (g.globalAlpha = this.defaultHighlightAlpha, g.fillStyle = this.defaultHighlightFill, g.fillRect((L.min - p) * this.scale, 0, (L.max - L.min) * this.scale, q));
                g.restore()
            }
            p = -1;
            
            //"center" == e.ruler ? p = 200 + (this.viewEnd - this.viewStart + 1) * this.scale / 2 : "left" == e.ruler ? p = 200 : "right" == e.ruler && (p = 200 + (this.viewEnd - this.viewStart + 1) * this.scale);
            p = 200 + (this.viewEnd - this.viewStart + 1) * this.scale / 2;
            
            0 <= p && (g.strokeStyle = "blue", g.beginPath(), g.moveTo(p, 0), g.lineTo(p, q), g.stroke());
            return z.toDataURL("image/png")
        }
    }, {
        "./cbrowser": 6,
        "./glyphs": 21,
        "./numformats": 26,
        "./sequence-draw": 31
    }],
    16: [function(e, u, s) {
        if ("undefined" !== typeof e) {
            var p = e("./cbrowser").Browser;
            e = e("./utils");
            var n = e.makeElement,
                m = e.removeChildren
        }
        p.prototype.openExportPanel =
            function() {
                var h = this;
                if ("export" === this.uiMode) this.hideToolPanel(), this.setUiMode("none");
                else {
                    var e = n("div", null, {
                            className: "export-form"
                        }),
                        p = n("select");
                    p.appendChild(n("option", "SVG", {
                        value: "svg"
                    }));
                    p.appendChild(n("option", "PNG", {
                        value: "png"
                    }));
                    /*
                    p.appendChild(n("option", "Dalliance config", {
                        value: "config"
                    }));
                    p.appendChild(n("option", "Dalliance sources", {
                        value: "sources"
                    }));
                    p.appendChild(n("option", "Dalliance page", {
                        value: "page"
                    }));
                    */
                    p.value = "svg";
                    p.addEventListener("change", function(a) {
                        m(g);
                        D()
                    }, !1);
                    e.appendChild(n("p", ["Export as: ", p]));
                    var z = n("input", null, {
                        type: "checkbox",
                        checked: this.exportHighlights
                    });
                    z.addEventListener("change", function(a) {
                        h.exportHighlights = z.checked;
                        h.storeStatus()
                    }, !1);
                    var a = n("input", null, {
                        type: "checkbox",
                        checked: this.exportRuler
                    });
                    a.addEventListener("change", function(b) {
                        h.exportRuler = a.checked;
                        h.storeStatus()
                    }, !1);
                    var b = n("input", null, {
                            type: "text",
                            value: "1.0"
                        }),
                        q = n("button", "Export", {
                            className: "btn btn-primary"
                        });
                    q.addEventListener("click", function(d) {
                        m(g);
                        var f,
                            l, q, e;
                        if ("svg" === p.value) f = URL.createObjectURL(h.makeSVG({
                            highlights: z.checked,
                            ruler: a.checked ? h.rulerLocation : "none"
                        })), l = "SVG", q = "image/svg", e = "dalliance-view.svg";
                        else if ("png" === p.value) {
                            q = parseFloat(b.value);
                            if (0.1 > q || 10 < q) {
                                alert("bad scale " + q);
                                return
                            }
                            f = h.exportImage({
                                highlights: z.checked,
                                ruler: a.checked ? h.rulerLocation : "none",
                                resolutionMultiplier: q
                            });
                            l = "Image";
                            q = "image/png";
                            e = "dalliance-view.png"
                        } else "config" === p.value ? (q = JSON.stringify(h.exportFullConfig(), null, 2), l = new Blob([q], {
                                type: "text/plain"
                            }),
                            f = URL.createObjectURL(l), l = "Configuration", q = "text/plain", e = "dalliance-config.json") : "sources" === p.value ? (q = JSON.stringify(h.exportSourceConfig(), null, 2), l = new Blob([q], {
                            type: "text/plain"
                        }), f = URL.createObjectURL(l), l = "Source array", q = "text/plain", e = "dalliance-sources.json") : "page" === p.value && (l = h.exportPageTemplate(), q = "text/html", l = new Blob([l], {
                            type: q
                        }), f = URL.createObjectURL(l), l = "Page template", e = "dalliance-view.html");
                        f && (d = n("a", "[Download]", {
                            href: f,
                            download: e,
                            type: q
                        }), q = n("a", "[Preview in browser]", {
                            href: f,
                            type: q,
                            target: "_new"
                        }), g.appendChild(n("p", ["" + l + " created: ", d, q])))
                    }, !1);
                    h.addViewListener(function() {
                        m(g)
                    });
                    h.addTierListener(function() {
                        m(g)
                    });
                    var g = n("p", ""),
                        l = n("tr", [n("th", "Include highlights", {}, {
                            width: "200px",
                            textAlign: "right"
                        }), n("td", z)]),
                        f = n("tr", [n("th", "Include vertical guideline"), n("td", a)]),
                        d = n("tr", [n("th", "Scale multiplier"), n("td", b)]),
                        t = n("table", [l, f, d]),
                        D = function() {
                            var a = p.value;
                            l.style.display = "svg" == a || "png" == a ? "table-row" : "none";
                            f.style.display = "svg" == a || "png" ==
                                a ? "table-row" : "none";
                            d.style.display = "png" == a ? "table-row" : "none"
                        };
                    D();
                    e.appendChild(t);
                    e.appendChild(q);
                    e.appendChild(g);
                    "none" !== this.uiMode && this.hideToolPanel();
                    this.browserHolder.insertBefore(e, this.svgHolder);
                    this.activeToolPanel = e;
                    this.setUiMode("export")
                }
            }
    }, {
        "./cbrowser": 6,
        "./utils": 48
    }],
    17: [function(e, u, s) {
        u = e("./cbrowser");
        s = e("./chainset");
        var p = e("./sourceadapters"),
            n = e("./utils");
        e = e("./das");
        window.Browser = u.Browser;
        window.sourcesAreEqual = u.sourcesAreEqual;
        window.Chainset = s.Chainset;
        window.makeElement = n.makeElement;
        window.dalliance_registerSourceAdapterFactory = p.registerSourceAdapterFactory;
        window.dalliance_registerParserFactory = p.registerParserFactory;
        window.dalliance_makeParser = p.makeParser;
        window.DASSequence = e.DASSequence;
        window.DASFeature = e.DASFeature;
        window.DASGroup = e.DASGroup;
        window.DASStylesheet = e.DASStylesheet;
        window.DASStyle = e.DASStyle;
        window.DASSource = e.DASSource
    }, {
        "./cbrowser": 6,
        "./chainset": 7,
        "./das": 10,
        "./sourceadapters": 34,
        "./utils": 48
    }],
    18: [function(e, u, s) {
        function p() {
            this.glyphs = [];
            this.height = 0;
            this.quant = null
        }

        function n(c) {
            Date.now();
            X = c.viewport.getContext("2d");
            b(c);
            c.padding = "number" === typeof c.dasSource.padding ? c.dasSource.padding : ca;
            var d = [],
                k = {},
                f = {},
                g;
            for (g in c.ungroupedFeatures)
                for (var l = c.ungroupedFeatures[g], q = 0; q < l.length; ++q) {
                    var e = l[q];
                    if (!e.parts) {
                        var w = c.styleForFeature(e);
                        if (w)
                            if ("LINEPLOT" == w.glyph) a(k, w.id, e), f[w.id] = w;
                            else {
                                var t = h(e, 0, w, c);
                                t && d.push(t)
                            }
                    }
                }
            for (var n in k) g = k[n], w = f[n], "LINEPLOT" == w.glyph && d.push(v(g, w, c));
            if (c.dasSource.collapseSuperGroups &&
                !c.bumped)
                for (var E in c.superGroups) {
                    k = c.superGroups[E];
                    c.groups[E] = z(c.groups[E]);
                    c.groups[E].isSuperGroup = !0;
                    f = {};
                    n = 1E10;
                    l = -1E10;
                    q = null;
                    for (t = 0; t < k.length; ++t)
                        if (g = c.groupedFeatures[k[t]]) {
                            for (w = 0; w < g.length; ++w) e = g[w], a(f, e.type, e), n = Math.min(e.min, n), l = Math.max(e.max, l), e.segment && !q && (q = e.segment);
                            if (c.groups[E] && !c.groups[E].links || 0 == c.groups[E].links.length) c.groups[E].links = c.groups[k[0]].links;
                            delete c.groupedFeatures[k[t]]
                        }
                    c.groups[E].max = l;
                    c.groups[E].min = n;
                    c.groups[E].segment = q;
                    for (var B in f) {
                        g =
                            f[B];
                        n = g[0];
                        l = null;
                        for (w = 0; w < g.length; ++w) e = g[w], e = new S(e.min, e.max), l = l ? K(l, e) : e;
                        l = l.ranges();
                        for (q = 0; q < l.length; ++q) {
                            for (var x = l[q], r = ((x.max() | 0) - (x.min() | 0) + 1) * k.length, y = 0, w = 0; w < g.length; ++w)
                                if (e = g[w], (e.min | 0) <= x.max() && (e.max | 0) >= x.min()) var D = Math.max(e.min | 0, x.min()),
                                    e = Math.min(e.max | 0, x.max()),
                                    y = y + (e - D + 1);
                            var e = new R,
                                A;
                            for (A in n) e[A] = n[A];
                            e.min = x.min();
                            e.max = x.max();
                            e.label && 1 < k.length && (e.label += " (" + k.length + " vars)");
                            e.visualWeight = 1 * y / r;
                            a(c.groupedFeatures, E, e)
                        }
                    }
                    delete c.superGroups[E]
                }
            A = [];
            for (var I in c.groupedFeatures) A.push(I);
            A.sort(function(a, b) {
                var d = c.groupedFeatures[a][0].score - c.groupedFeatures[b][0].score;
                return 0 < d ? -1 : 0 == d ? 0 : 1
            });
            B = {};
            for (e = 0; e < A.length; ++e)
                if (I = A[e], t = m(c.groupedFeatures[I], 0, c.groups[I], c, c.dasSource.collapseSuperGroups && !c.bumped ? "collapsed_gene" : "tent")) t.group = c.groups[I], B[I] = t;
            for (E in c.superGroups) {
                k = c.superGroups[E];
                I = [];
                n = 1E10;
                l = -1E10;
                for (A = 0; A < k.length; ++A) e = B[k[A]], B[k[A]] = null, e && (I.push(e), n = Math.min(n, e.min()), l = Math.max(l, e.max()));
                for (A =
                    0; A < I.length; ++A) e = I[A], d.push(new G(e, n, l))
            }
            for (t in B)(e = B[t]) && d.push(e);
            E = new p;
            I = [];
            B = c.subtierMax || c.dasSource.subtierMax || c.browser.defaultSubtierMax;
            A = !1;
            e = 0;
            a: for (; e < d.length; ++e)
                if (t = d[e], t.bump && (c.bumped || c.dasSource.collapseSuperGroups)) {
                    for (k = 0; k < I.length; ++k)
                        if (f = I[k], f.hasSpaceFor(t)) {
                            f.add(t);
                            continue a
                        }
                    I.length >= B ? A = !0 : (f = new p, f.add(t), I.push(f))
                } else E.add(t);
            0 < E.glyphs.length && (I = [E].concat(I));
            for (k = 0; k < I.length; ++k) f = I[k], f.quant && f.glyphs.unshift(new C(f.height));
            for (k = 0; k <
                I.length; ++k) f = I[k], f.glyphs.sort(function(a, c) {
                return (a.zindex || 0) - (c.zindex || 0)
            });
            c.subtiers = I;
            c.glyphCacheOrigin = c.browser.viewStart;
            A ? c.updateStatus("Bumping limit exceeded, use the track editor to see more features") : c.updateStatus()
        }

        function m(a, c, b, d, k) {
            k = d.styleForFeature(b);
            var l, e = !1,
                q = [];
            c = null;
            for (var w = 0; w < a.length; ++w) {
                var t = a[w];
                t.orientation && null == c && (c = t.orientation);
                !l && t.label && (l = t.label);
                var m = d.styleForFeature(t);
                m && !t.parts && (V(m.LABEL) && (e = !0), (t = h(t, 0, m, d, null, !0)) && q.push(t))
            }
            if (0 ==
                q.length) return null;
            a = "flat";
            k && "LINE" === k.glyph || (d.dasSource.collapseSuperGroups && !d.bumped ? "+" === c ? a = "collapsed+" : "-" === c && (a = "collapsed-") : "+" === c ? a = "hat+" : "-" === c && (a = "hat-"));
            d = null;
            if (l && e || k && (V(k.LABEL) || V(k.LABELS))) d = b.label || l;
            b = new g(q, a);
            d && ("+" === c ? d = ">" + d : "-" === c && (d = "<" + d), b = new f(X, b, d, !1));
            b.bump = !0;
            return b
        }

        function h(a, b, k, l, h, e) {
            function m(a, c, b) {
                var k = null;
                if (a.currentSequence) {
                    var d = a.currentSequence.start | 0,
                        f = a.currentSequence.end | 0;
                    if (d <= b && f >= c) {
                        for (var g = Math.max(c, d),
                                f = Math.min(b, f), k = a.currentSequence.seq.substr(g - d, f - g + 1); c < g;) k = "N" + k, g--;
                        for (; b > f;) k += "N", f++
                    }
                }
                return k
            }
            var n = l.browser.scale,
                p = l.browser.viewStart,
                r = k.glyph || "BOX",
                z = a.min,
                C = a.max,
                s = a.orientation,
                K = a.score,
                M = a.label || a.id,
                v = (z - p) * n,
                S = Math.max((C - p + 1) * n, v + 1),
                R = l.forceHeight || k.HEIGHT || h || 12,
                u = R *= 1,
                p = k.BUMP && V(k.BUMP),
                T;
            if ("CROSS" === r || "EX" === r || "TRIANGLE" === r || "DOT" === r || "SQUARE" === r || "STAR" === r || "PLIMSOLL" === r) {
                var W = k.FGCOLOR || "black",
                    O = k.BGCOLOR || "none",
                    z = k.STROKECOLOR;
                k.BGITEM && a.itemRgb ? W =
                    a.itemRgb : V(k.COLOR_BY_SCORE2) && (n = k.BGGRAD || k._gradient, n || (n = E(50, k.COLOR1, k.COLOR2, k.COLOR3), k._gradient = n), R = a.score2, void 0 == R && W || (u = k.MIN2 ? 1 * k.MIN2 : 0, C = (1 * (R || 0) - u) / ((k.MAX2 ? 1 * k.MAX2 : 1) - u) * n.length | 0, 0 > C && (C = 0), C >= n.length && (C = n.length - 1), W = n[C]));
                R = l.forceHeight || k.HEIGHT || h || 12;
                u = R *= 1;
                C = k.SIZE || R;
                k.RSIZE && (C = 1 * k.RSIZE * R);
                k.STROKETHRESHOLD && C < 1 * k.STROKETHRESHOLD && (z = null);
                C *= 1;
                R = (v + S) / 2;
                n = C / 2;
                "EX" === r ? z = new t(R, C, W) : "TRIANGLE" === r ? z = new D(R, C, k.DIRECTION || "N", k.LINEWIDTH || C, W, z) : "DOT" === r ?
                    z = new P(R, C, W, z) : "PLIMSOLL" === r ? z = new J(R, C, 0.2 * C, W, z) : "SQUARE" === r ? z = new q(R - n, 0, C, C, W, z) : "STAR" === r ? (r = 5, k.POINTS && (r = k.POINTS | 0), z = new A(R, n, r, W, z)) : z = new d(R, C, W);
                O && "none" != O && 5 < S - v && (v = new q(v, 0, S - v, C, O), z = new g([v, z]));
                if (V(k.SCATTER)) {
                    O = l.quantMin(k);
                    (C = l.quantMax(k)) || (C = 0 > O ? 0 : 10);
                    O || (O = 0);
                    l = (1 * K - O) / (C - O);
                    K = -1 * O / (C - O);
                    if (0 > l || 1 < l) return null;
                    l >= K ? (R = Math.max(1, (l - K) * u), b = b + (1 - K) * u - R) : (R = Math.max(1, (l - K) * u), b += (1 - K) * u);
                    T = {
                        min: O,
                        max: C
                    };
                    v = 0;
                    S = "undefined" !== typeof a.forceLabel ? a.forceLabel : k.LABEL;
                    Z(S) && M && !e && (z = new f(X, z, M, !0, null, "above" == S ? "above" : "below"), "above" == S && (v = z.textHeight + 2), e = !0);
                    z = new I(z, 0, b - n - v, u)
                }
            } else if ("HISTOGRAM" === r || "GRADIENT" === r && "undefined" !== K) O = l.quantMin(k), (C = l.quantMax(k)) || (C = 0 > O ? 0 : 10), O || (O = 0), 1 * K < 1 * O && (K = O), 1 * K > 1 * C && (K = C), l = (1 * K - O) / (C - O), K = -1 * O / (C - O), "HISTOGRAM" === r && (l >= K ? (R = (l - Math.max(0, K)) * u, b = b + (1 - Math.max(0, K)) * u - R) : (R = (Math.max(0, K) - l) * u, b += (1 - Math.max(0, K)) * u), T = {
                    min: O,
                    max: C
                }), W = k.FGCOLOR || null, O = k.BGCOLOR || k.COLOR1 || "green", k.BGITEM && a.itemRgb &&
                (O = a.itemRgb), K = k.ALPHA ? 1 * k.ALPHA : null, k.BGGRAD && (n = k.BGGRAD, C = l * n.length | 0, 0 > C && (C = 0), C >= n.length && (C = n.length - 1), O = n[C]), k.COLOR2 && (n = k._gradient, n || (n = E(50, k.COLOR1, k.COLOR2, k.COLOR3), k._gradient = n), C = l * n.length | 0, 0 > C && (C = 0), C >= n.length && (C = n.length - 1), O = n[C]), z = new q(v, b, S - v, R, O, W, K), z = new I(z, 0, 0, u);
            else if ("HIDDEN" === r) z = new G(null, v, S), e = !0;
            else if ("ARROW" === r) l = k.FGCOLOR || "purple", b = V(k.PARALLEL), u = V(k.SOUTHWEST), K = V(k.NORTHEAST), z = new H(v, S, R, l, b, u, K);
            else if ("ANCHORED_ARROW" === r) W = k.FGCOLOR ||
                "none", O = k.BGCOLOR || "green", z = new L(v, S, R, O, W, s), z.bump = !0;
            else if ("SPAN" === r) W = k.FGCOLOR || "black", z = new N(v, S, R, W);
            else if ("LINE" === r) W = k.FGCOLOR || "black", z = new Q(v, S, R, k.STYLE || "solid", s, W);
            else if ("PRIMERS" === r) W = k.FGCOLOR || "black", O = k.BGCOLOR || "red", z = new y(v, S, R, O, W);
            else if ("TEXT" === r) l = k.STRING || "text", O = k.FGCOLOR || "black", z = new c(X, v, S, R, O, l);
            else if ("TOOMANY" === r) W = k.FGCOLOR || "gray", O = k.BGCOLOR || "orange", z = new F(v, S, R, O, W);
            else if ("POINT" === r) R = l.forceHeight || k.HEIGHT || 30, O = l.quantMin(k),
                C = l.quantMax(k), l = (1 * K - O) / (C - O), b = 1 * R / (C - O) * (K - 1 * O) | 0, T = {
                    min: O,
                    max: C
                }, O = k.FGCOLOR || k.COLOR1 || "black", k.COLOR2 && (n = k._gradient, n || (n = E(50, k.COLOR1, k.COLOR2, k.COLOR3), k._gradient = n), C = l * n.length | 0, 0 > C && (C = 0), C >= n.length && (C = n.length - 1), O = n[C]), z = new x((v + S) / 2, R - b, R, O);
            else if ("__SEQUENCE" === r) {
                K = r = a.seq;
                O = W = a.quals;
                b = V(k.__INSERTIONS);
                u = [];
                if (a.cigar) {
                    h = ja(a.cigar);
                    for (var O = K = "", ca = s = 0; ca < h.length; ++ca) {
                        var aa = h[ca];
                        if ("M" == aa.op) K += r.substr(s, aa.cnt), O += W.substr(s, aa.cnt), s += aa.cnt;
                        else if ("D" == aa.op)
                            for (var ia =
                                    0; ia < aa.cnt; ++ia) K += "-", O += "Z";
                        else if ("I" == aa.op) {
                            var ba = r.substr(s, aa.cnt),
                                ia = new D(v + K.length * n, 5, "S", 5, l.browser.baseColors.I);
                            b && (ia = new f(X, ia, ba, !1, "center", "above", "7px sans-serif"));
                            ia.feature = {
                                label: "Insertion: " + ba,
                                type: "insertion",
                                method: "insertion"
                            };
                            u.push(ia);
                            s += aa.cnt
                        } else "S" == aa.op ? s += aa.cnt : console.log("unknown cigop" + aa.op)
                    }
                }
                n = m(l, z, C);
                if (K && n && ("mismatch" === k.__SEQCOLOR || "mismatch-all" === k.__SEQCOLOR)) {
                    z = [];
                    C = "-" === a.orientation ? "," : ".";
                    for (r = 0; r < K.length; ++r) z.push(K[r] == n[r] ?
                        C : K[r]);
                    K = z.join("")
                }
                z = "-" === a.orientation ? k._minusColor || "lightskyblue" : k._plusColor || "lightsalmon";
                k.__disableQuals && (O = !1);
                z = new w(l.browser.baseColors, z, v, S, R, K, n, k.__SEQCOLOR, O, !V(k.__CLEARBG));
                b && (z = new I(z, 0, 7));
                0 < u.length && (u.splice(0, 0, z), z = new g(u))
            } else if ("__INSERTION" === r) ia = new D(v, 5, "S", 5, l.browser.baseColors.I), z = new f(X, ia, a.insertion || a.altAlleles[0], !1, "center", "above", "7px sans-serif"), 1 < S - v && (O = k.BGCOLOR || k.COLOR1 || "green", v = new q(v, 5, S - v, R, O, W), z = new g([v, z]));
            else {
                if ("__NONE" ===
                    r) return null;
                W = k.FGCOLOR || null;
                O = k.BGCOLOR || k.COLOR1 || "green";
                k.BGITEM && a.itemRgb && (O = a.itemRgb);
                n = (S - v) / (C - z);
                "translation" == a.type && ("protein_coding" == a.method || a.readframeExplicit) && (!a.tags || 0 > a.tags.indexOf("cds_start_NF") || a.readframeExplicit) && (!l.dasSource.collapseSuperGroups || l.bumped) && 0.5 <= n ? (n = m(l, z, C), z = new B(v, S, R, O, n, a.orientation, a.readframe)) : z = new q(v, 0, S - v, R, O, W)
            }(V(k.LABEL) || a.forceLabel) && M && !e && (z = new f(X, z, M, !1));
            p && (z.bump = !0);
            z.feature = a;
            T && (z.quant = T);
            k.ZINDEX && (z.zindex =
                k.ZINDEX | 0);
            return z
        }

        function v(a, c, b) {
            var k = b.browser.viewStart,
                d = b.browser.scale,
                f = b.forceHeight || c.HEIGHT || 30,
                g = b.quantMin(c);
            b = b.quantMax(c);
            for (var h = 1 * f / (b - g), e = c.FGCOLOR || c.COLOR1 || "black", q = [], w = 0; w < a.length; ++w) {
                var t = a[w],
                    n = f - ((t.score - 1 * g) * h | 0);
                q.push((((t.min | 0) + (t.max | 0)) / 2 - k) * d);
                q.push(n)
            }
            a = new l(q, e, f);
            a.quant = {
                min: g,
                max: b
            };
            c.ZINDEX && (a.zindex = c.ZINDEX | 0);
            return a
        }
        if ("undefined" !== typeof e) {
            var r = e("./utils"),
                z = r.shallowCopy,
                a = r.pusho,
                r = e("./tier").DasTier,
                b = e("./features").sortFeatures;
            s = e("./glyphs");
            var q = s.BoxGlyph,
                g = s.GroupGlyph,
                l = s.LineGraphGlyph,
                f = s.LabelledGlyph,
                d = s.CrossGlyph,
                t = s.ExGlyph,
                D = s.TriangleGlyph,
                P = s.DotGlyph,
                G = s.PaddedGlyph,
                L = s.AArrowGlyph,
                N = s.SpanGlyph,
                Q = s.LineGlyph,
                y = s.PrimersGlyph,
                H = s.ArrowGlyph,
                F = s.TooManyGlyph,
                c = s.TextGlyph,
                w = s.SequenceGlyph,
                B = s.AminoAcidGlyph,
                I = s.TranslatedGlyph,
                x = s.PointGlyph,
                C = s.GridGlyph,
                A = s.StarGlyph,
                J = s.PlimsollGlyph,
                k = s.OverlayLabelCanvas,
                E = e("./color").makeGradient;
            s = e("./spans");
            var S = s.Range,
                K = s.union;
            s = e("./das");
            var R = s.DASFeature,
                V = s.isDasBooleanTrue,
                Z = s.isDasBooleanNotFalse,
                ja = e("./cigar").parseCigar,
                T = e("./numformats").formatQuantLabel
        }
        var ca = 3;
        p.prototype.indexFor = function(a) {
            a = a.min();
            for (var c = 0, b = this.glyphs.length; b > c;) {
                var k = (c + b) / 2 | 0;
                if (k >= this.glyphs.length) return this.glyphs.length;
                a < this.glyphs[k].min() ? b = k : c = k + 1
            }
            return b
        };
        p.prototype.add = function(a) {
            var c = this.indexFor(a);
            this.glyphs.splice(c, 0, a);
            this.height = Math.max(this.height, a.height());
            a.quant && null == this.quant && (this.quant = a.quant)
        };
        p.prototype.hasSpaceFor =
            function(a) {
                var c = this.indexFor(a);
                return 0 < c && this.glyphs[c - 1].max() >= a.min() || c < this.glyphs.length && this.glyphs[c].min() <= a.max() ? !1 : !0
            };
        var X;
        r.prototype.paint = function() {
            var a = this.browser.retina && (1 < window.devicePixelRatio || 1 < this.browser.devicePixelRatio),
                c = this.subtiers;
            if (c) {
                var b = this.browser.featurePanelWidth + 2E3;
                a && (b *= 2);
                var d = this.viewport.width | 0;
                d < b - 50 && (this.viewport.width = d = b);
                for (var b = this.padding, f = 0; f < c.length; ++f) b = b + c[f].height + this.padding;
                f = b = Math.max(b + 6, this.browser.minTierHeight);
                a && (f *= 2);
                f != this.viewport.height && (this.viewport.height = f);
                Math.max(b, this.browser.minTierHeight);
                this.viewportHolder.style.left = "-1000px";
                this.viewport.style.width = a ? "" + d / 2 + "px" : "" + d + "px";
                this.viewport.style.height = "" + b + "px";
                this.layoutHeight = Math.max(b, this.browser.minTierHeight);
                this.updateHeight();
                this.norigin = this.browser.viewStart;
                c = this.viewport.getContext("2d");
                c.clearRect(0, 0, d, f);
                c.save();
                a && c.scale(2, 2);
                d = this.browser.viewStart - 1E3 / this.browser.scale;
                f = this.browser.viewEnd + 1E3 / this.browser.scale;
                a = [];
                if (this.knownCoverage)
                    for (var g = this.knownCoverage.ranges(), l = 0; l < g.length; ++l) {
                        var h = g[l];
                        0 == l ? h.min() > d && a.push({
                            min: d,
                            max: h.min() - 1
                        }) : a.push({
                            min: g[l - 1].max() + 1,
                            max: h.min() - 1
                        });
                        l == g.length - 1 && h.max() < f && a.push({
                            min: h.max() + 1,
                            max: f
                        })
                    }
                if (0 < a.length)
                    for (c.fillStyle = "gray", d = 0; d < a.length; ++d) f = a[d], g = (f.min - this.browser.viewStart) * this.browser.scale + 1E3, c.fillRect(g, 0, (f.max - this.browser.viewStart) * this.browser.scale + 1E3 - g, b);
                b = new k;
                a = (this.glyphCacheOrigin - this.browser.viewStart) * this.browser.scale +
                    1E3;
                c.translate(a, this.padding);
                b.translate(0, this.padding);
                this.paintToContext(c, b, a);
                this.overlayLabelCanvas = 0 < b.glyphs.length ? b : null;
                c.restore();
                this.drawOverlay();
                this.paintQuant()
            }
        };
        r.prototype.paintToContext = function(a, c, b) {
            var k = this.subtiers,
                d = this.viewport.width | 0;
            a.save();
            for (var f = 0; f < k.length; ++f) {
                for (var g = null, l = k[f].glyphs, h = 0; h < l.length; ++h) {
                    var e = l[h];
                    e.min() < d - b && e.max() > -b && (e = l[h], e.draw(a, c), e.quant && (g = e.quant))
                }
                a.translate(0, k[f].height + this.padding);
                c.translate(0, k[f].height +
                    this.padding)
            }
            a.restore();
            g && this.quantLeapThreshold && this.featureSource && this.browser.sourceAdapterIsCapable(this.featureSource, "quantLeap") && (c = k[0].height * (1 - (this.quantLeapThreshold - g.min) / (g.max - g.min)), a.save(), a.strokeStyle = "red", a.lineWidth = 0.3, a.beginPath(), a.moveTo(-1E3, c), a.lineTo(d + 1E3, c), a.stroke(), a.restore())
        };
        r.prototype.paintQuant = function() {
            if (this.quantOverlay) {
                var a = this.browser.retina && (1 < window.devicePixelRatio || 1 < this.browser.devicePixelRatio),
                    c;
                this.subtiers && 0 < this.subtiers.length &&
                    (c = this.subtiers[0].quant);
                if (c) {
                    var b = this.subtiers[0].height;
                    this.quantOverlay.height = this.viewport.height;
                    this.quantOverlay.width = a ? 100 : 50;
                    this.quantOverlay.style.height = "" + (a ? this.quantOverlay.height / 2 : this.quantOverlay.height) + "px";
                    this.quantOverlay.style.width = "50px";
                    this.quantOverlay.style.display = "block";
                    var k = this.quantOverlay.getContext("2d");
                    a && k.scale(2, 2);
                    a = 2;
                    40 < b && (a = 1 + (b / 20 | 0));
                    var d = (b + 2 * this.padding) / (a - 1),
                        f = (c.max - c.min) / (a - 1);
                    k.fillStyle = "white";
                    k.globalAlpha = 0.6;
                    "right" == this.browser.rulerLocation ?
                        k.fillRect(20, 0, 30, b + 2 * this.padding) : k.fillRect(0, 0, 30, b + 2 * this.padding);
                    k.globalAlpha = 1;
                    k.strokeStyle = "black";
                    k.lineWidth = 1;
                    k.beginPath();
                    if ("right" == this.browser.rulerLocation) {
                        k.moveTo(42, this.padding);
                        k.lineTo(50, this.padding);
                        k.lineTo(50, b + this.padding);
                        k.lineTo(42, b + this.padding);
                        for (var g = 1; g < a - 1; ++g) {
                            var l = g * d;
                            k.moveTo(50, l);
                            k.lineTo(45, l)
                        }
                    } else
                        for (k.moveTo(8, this.padding), k.lineTo(0, this.padding), k.lineTo(0, b + this.padding), k.lineTo(8, b + this.padding), g = 1; g < a - 1; ++g) l = g * d, k.moveTo(0, l), k.lineTo(5,
                            l);
                    k.stroke();
                    k.fillStyle = "black";
                    if ("right" == this.browser.rulerLocation)
                        for (k.textAlign = "right", k.fillText(T(c.max), 41, 8), k.fillText(T(c.min), 41, b + this.padding), g = 1; g < a - 1; ++g) l = g * d, k.fillText(T(1 * c.max - g * f), 41, l + 3);
                    else
                        for (k.textAlign = "left", k.fillText(T(c.max), 9, 8), k.fillText(T(c.min), 9, b + this.padding), g = 1; g < a - 1; ++g) l = g * d, k.fillText(T(1 * c.max - g * f), 9, l + 3)
                } else this.quantOverlay.style.display = "none"
            }
        };
        r.prototype.styleForFeature = function(a) {
            var c = this.browser.zoomForCurrentScale();
            if (!this.stylesheet) return null;
            for (var b = null, k = this.stylesheet.styles, d = 0; d < k.length; ++d) {
                var f = k[d];
                if (!(f.zoom && f.zoom != c || f.orientation && f.orientation != a.orientation)) {
                    var g = f._labelRE;
                    g && g.test || (g = new RegExp("^" + f.label + "$"), f._labelRE = g);
                    if (!f.label || g.test(a.label))
                        if (g = f._methodRE, g && g.test || (g = new RegExp("^" + f.method + "$"), f._methodRE = g), !f.method || g.test(a.method)) {
                            if (f.type)
                                if ("default" == f.type) {
                                    b || (b = f.style);
                                    continue
                                } else if (g = f._typeRE, g && g.test || (g = new RegExp("^" + f.type + "$"), f._typeRE = g), !g.test(a.type)) continue;
                            return f.style
                        }
                }
            }
            return b
        };
        r.prototype.quantMin = function(a) {
            return this.forceMinDynamic ? this.currentFeaturesMinScore || 0 : "number" === typeof this.forceMin ? this.forceMin : a.MIN || this.currentFeaturesMinScore || 0
        };
        r.prototype.quantMax = function(a) {
            return this.forceMaxDynamic ? this.currentFeaturesMaxScore || 0 : "number" === typeof this.forceMax ? this.forceMax : a.MAX || this.currentFeaturesMaxScore || 0
        };
        "undefined" !== typeof u && (u.exports = {
            drawFeatureTier: n
        })
    }, {
        "./cigar": 8,
        "./color": 9,
        "./das": 10,
        "./features": 20,
        "./glyphs": 21,
        "./numformats": 26,
        "./spans": 35,
        "./tier": 44,
        "./utils": 48
    }],
    19: [function(e, u, s) {
        function p(a, b, e) {
            var g = h(e.type, b.type),
                l = h(e.label, b.label, e.id, b.id);
            l && 0 != l.indexOf("__dazzle") && (g = g + ": " + l);
            this.hit = a;
            this.feature = b;
            this.group = e;
            this.title = g;
            this.sections = []
        }

        function n(a, b) {
            var h = [];
            if (a)
                for (var g = 0; g < a.length; ++g) v(h, a[g]);
            if (b)
                for (g = 0; g < b.length; ++g) v(h, b[g]);
            return h
        }
        if ("undefined" !== typeof e) {
            var m = e("./cbrowser").Browser;
            e = e("./utils");
            var h = e.pick,
                v = e.pushnew,
                r = e.makeElement
        }
        var z = /^([A-Za-z_-]+)=(.+)/;
        m.prototype.addFeatureInfoPlugin = function(a) {
            this.featureInfoPlugins || (this.featureInfoPlugins = []);
            this.featureInfoPlugins.push(a)
        };
        p.prototype.setTitle = function(a) {
            this.title = a
        };
        p.prototype.add = function(a, b) {
            "string" === typeof b && (b = r("span", b));
            this.sections.push({
                label: a,
                info: b
            })
        };
        m.prototype.featurePopup = function(a, b, h, g) {
            var l = h.length;
            b = 0 <= --l ? h[l] : {};
            var f = 0 <= --l ? h[l] : {};
            h = new p(h, b, f);
            h.tier = g;
            for (var l = this.featureInfoPlugins || [], d = 0; d < l.length; ++d) try {
                l[d](b, h)
            } catch (e) {
                console.log(e.stack ||
                    e)
            }
            l = g.featureInfoPlugins || [];
            for (d = 0; d < l.length; ++d) try {
                l[d](b, h)
            } catch (m) {
                console.log(m.stack || m)
            }
            this.removeAllPopups();
            g = r("table", null, {
                className: "table table-striped table-condensed"
            });
            g.style.width = "100%";
            g.style.margin = "0px";
            l = 0;
            b.method && (d = r("tr", [r("th", "Method"), r("td", b.method)]), g.appendChild(d), ++l);
            d = f.segment ? f : b;
            d = r("tr", [r("th", "Location"), r("td", d.segment + ":" + d.min + "-" + d.max, {}, {
                minWidth: "200px"
            })]);
            g.appendChild(d);
            ++l;
            void 0 === b.score || null === b.score || "-" == b.score || b.suppressScore ||
                (d = r("tr", [r("th", "Score"), r("td", "" + b.score)]), g.appendChild(d), ++l);
            (d = n(f.links, b.links)) && 0 < d.length && (d = r("tr", [r("th", "Links"), r("td", d.map(function(a) {
                return r("div", r("a", a.desc, {
                    href: a.uri,
                    target: "_new"
                }))
            }))]), g.appendChild(d), ++l);
            b = n(f.notes, b.notes);
            for (f = 0; f < b.length; ++f) {
                var d = "Note",
                    v = b[f],
                    G = v.match(z);
                G && (d = G[1], v = G[2]);
                d = r("tr", [r("th", d), r("td", v)]);
                g.appendChild(d);
                ++l
            }
            for (b = 0; b < h.sections.length; ++b) l = h.sections[b], g.appendChild(r("tr", [r("th", l.label), r("td", l.info)]));
            this.popit(a,
                h.title || "Feature", g, {
                    width: 450
                })
        }
    }, {
        "./cbrowser": 6,
        "./utils": 48
    }],
    20: [function(e, u, s) {
        function p(h) {
            for (var e = h.browser.drawnStart, p = h.browser.drawnEnd, z = {}, a = {}, b = {}, q = {}, g = {}, l = {}, f = {}, d = {}, t = [], D, P, G, L = function() {
                    G = {};
                    for (var a = 0; a < h.currentFeatures.length; ++a) {
                        var c = h.currentFeatures[a];
                        c.id && (G[c.id] = c)
                    }
                }, N = function(a) {
                    var c = [];
                    if (a.parents)
                        for (var b = 0; b < a.parents.length; ++b) {
                            var d = a.parents[b],
                                f = G[d];
                            f && "SO:0000704" == f.typeCv && pushnew(c, d)
                        }
                    return c
                }, Q = 0; Q < h.currentFeatures.length; ++Q) {
                var y =
                    h.currentFeatures[Q];
                if (!y.parts) {
                    var H = y.min <= p && y.max >= e;
                    if (y.min && y.max) {
                        if (y.score && "." != y.score && "-" != y.score) {
                            var s = 1 * y.score;
                            if (!D || s < D) D = s;
                            if (!P || s > P) P = s
                        }
                        var s = [],
                            c = null;
                        if (y.groups)
                            for (var w = 0; w < y.groups.length; ++w) {
                                var B = y.groups[w],
                                    I = B.id;
                                if ("gene" == B.type) c = I, l[I] = B;
                                else if ("translation" != B.type) {
                                    n(a, I, y);
                                    l[I] = B;
                                    s.push(I);
                                    B = q[I];
                                    if (!B || y.min < B) q[I] = y.min;
                                    B = g[I];
                                    if (!B || y.max > B) g[I] = y.max
                                }
                            }
                        if (y.parents)
                            for (G || L(), w = 0; w < y.parents.length; ++w) {
                                var x = y.parents[w],
                                    C = G[x];
                                if (C) {
                                    C.parts || (C.parts = [y]);
                                    m(a, x, C);
                                    n(a, x, y);
                                    l[x] || (l[x] = {
                                        type: C.type,
                                        id: C.id,
                                        label: C.label || C.id
                                    });
                                    s.push(x);
                                    B = q[x];
                                    if (!B || y.min < B) q[x] = y.min;
                                    B = g[x];
                                    if (!B || y.max > B) g[x] = y.max;
                                    B = N(C);
                                    0 < B.length && (c = B[0], x = G[B[0]], l[B[0]] = {
                                        type: x.type,
                                        id: x.id,
                                        label: x.label || x.id
                                    }, h.dasSource.collapseSuperGroups || (h.dasSource.collapseSuperGroups = !0))
                                }
                            }
                        if (0 == s.length) H && n(z, y.type, y);
                        else if (c)
                            for (B = 0; B < s.length; ++B) I = s[B], m(f, c, I), d[I] = c
                    } else t.push(y)
                }
            }
            for (I in a) t = l[I], "number" !== typeof t.min && (t.min = q[I]), "number" !== typeof t.max && (t.max =
                g[I]), g[I] >= e && q[I] <= p && (b[I] = a[I]);
            h.ungroupedFeatures = z;
            h.groupedFeatures = b;
            h.groups = l;
            h.superGroups = f;
            h.groupsToSupers = d;
            D && (0 < D ? D = 0 : 0 > P && (P = 0), h.currentFeaturesMinScore = D, h.currentFeaturesMaxScore = P)
        }
        if ("undefined" !== typeof e) {
            e = e("./utils");
            var n = e.pusho,
                m = e.pushnewo
        }
        "undefined" !== typeof u && (u.exports = {
            sortFeatures: p
        })
    }, {
        "./utils": 48
    }],
    21: [function(e, u, s) {
        function p(a, c) {
            this._stroke = a;
            this._fill = c
        }

        function n(a, c, b, d, f, g, l, h) {
            this.x = a;
            this.y = c;
            this._width = b;
            this._height = d;
            this.fill = f;
            this.stroke =
                g;
            this._alpha = l;
            this._radius = h || 0
        }

        function m(a, c) {
            this.glyphs = a;
            this.connector = c;
            this.h = a[0].height();
            for (var b = [], d = 0; d < a.length; ++d) {
                var f = a[d];
                b.push(new I(f.min(), f.max()));
                this.h = Math.max(this.h, f.height())
            }
            this.coverage = B(b)
        }

        function h(a, c, b) {
            this.points = a;
            this.color = c;
            this._height = b || 50
        }

        function v(a, c, b, d, f, g, l) {
            this.glyph = c;
            this.text = b;
            this.anchor = f || "left";
            this.align = g || "below";
            l && (this.font = l);
            this.font && (a.save(), a.font = this.font);
            b = a.measureText(b);
            this.font && a.restore();
            this.textLen =
                b.width;
            this.textHeight = 5;
            this.bump = c.bump;
            this.measured = !d
        }

        function r(a, c, b) {
            this._x = a;
            this._height = c;
            this._stroke = b
        }

        function z(a, c, b) {
            this._x = a;
            this._height = c;
            this._stroke = b
        }

        function a(a, c, b, d, f, g) {
            p.call(this, g, f);
            this._x = a;
            this._height = c;
            this._dir = b;
            this._width = d
        }

        function b(a, c, b, d) {
            this._x = a;
            this._height = c;
            this._fill = b;
            this._stroke = d
        }

        function q(a, c, b) {
            this.glyph = a;
            this._min = c;
            this._max = b;
            a && (this.bump = a.bump)
        }

        function g(a, c, b, d, f, g) {
            p.call(this, f, d);
            this._min = a;
            this._max = c;
            this._height = b;
            this._ori =
                g
        }

        function l(a, c, b, d) {
            p.call(this, d, null);
            this._min = a;
            this._max = c;
            this._height = b
        }

        function f(a, c, b, d, f, g) {
            this._min = a;
            this._max = c;
            this._height = b;
            this._style = d;
            this._strand = f;
            this._stroke = g
        }

        function d(a, c, b, d, f) {
            this._min = a;
            this._max = c;
            this._height = b;
            this._fill = d;
            this._stroke = f
        }

        function t(a, c, b, d, f, g, l) {
            p.call(this, null, d);
            this._min = a;
            this._max = c;
            this._height = b;
            this._color = d;
            this._parallel = f;
            this._sw = g;
            this._ne = l
        }

        function D(a, c, b, d, f) {
            this._min = a;
            this._max = c;
            this._height = b;
            this._fill = d;
            this._stroke =
                f
        }

        function P(a, c, b, d, f, g) {
            this._min = c;
            this._max = b;
            this._height = d;
            this._fill = f;
            this._string = g;
            this._textLen = a.measureText(g).width
        }

        function G(a, c, b) {
            var d = {
                red: "darkred",
                purple: "mediumpurple",
                blue: "darkblue",
                green: "darkgreen"
            }[b.toLowerCase()];
            b = d ? [b, d] : ["rgb(73, 68, 149)", "rgb(9, 0, 103)"];
            return "?" == a ? "black" : "M" == a ? "greenyellow" : "*" == a ? "crimson" : b[c % 2]
        }

        function L(a) {
            var c = {
                A: "T",
                T: "A",
                G: "C",
                C: "G"
            };
            a = a.split("").reverse().join("");
            for (var b = [], d = 0; d < a.length; ++d) {
                var f = a.substr(d, 1).toUpperCase();
                b.push(f in c ? c[f] : "N")
            }
            return b.join("")
        }

        function N(a, c, b, d, f, g, l) {
            this._min = a;
            this._max = c;
            this._height = b;
            this._fill = d;
            this._seq = f;
            this._orientation = g;
            this._readframe = l
        }

        function Q(a, c, b, d) {
            this.glyph = a;
            this._height = d;
            this._x = c;
            this._y = b
        }

        function y(a, c, b, d) {
            this._x = a;
            this._y = c;
            this._height = b;
            this._fill = d
        }

        function H(a) {
            this._height = a || 50
        }

        function F(a, c, b, d, f) {
            p.call(this, f, d);
            this._x = a;
            this._r = c;
            this._points = b
        }

        function c(a, c, b, d, f) {
            this._x = a;
            this._height = c;
            this._overhang = b;
            this._fill = d;
            this._stroke =
                f;
            this._hh = c / 2
        }

        function w() {
            this.oy = this.ox = 0;
            this.glyphs = []
        }
        if ("undefined" !== typeof e) {
            s = e("./spans");
            var B = s.union,
                I = s.Range;
            s = e("./utils");
            var x = s.makeElementNS,
                C = s.AMINO_ACID_TRANSLATION;
            e = e("./svg-utils");
            var A = e.NS_SVG,
                J = e.SVGPath
        }
        p.prototype.draw = function(a) {
            a.beginPath();
            this.drawPath(a);
            this._fill && (a.fillStyle = this._fill, a.fill());
            this._stroke && (a.strokeStyle = this._stroke, a.stroke())
        };
        p.prototype.toSVG = function() {
            var a = new J;
            this.drawPath(a);
            return x(A, "path", null, {
                d: a.toPathData(),
                fill: this._fill ||
                    "none",
                stroke: this._stroke || "none"
            })
        };
        p.prototype.drawPath = function(a) {
            throw "drawPath method on PathGlyphBase must be overridden";
        };
        n.prototype.draw = function(a) {
            var c = this._radius;
            0 < c ? (a.beginPath(), a.moveTo(this.x + c, this.y), a.lineTo(this.x + this._width - c, this.y), a.arcTo(this.x + this._width, this.y, this.x + this._width, this.y + c, c), a.lineTo(this.x + this._width, this.y + this._height - c), a.arcTo(this.x + this._width, this.y + this._height, this.x + this._width - c, this.y + this._height, c), a.lineTo(this.x + c, this.y + this._height),
                a.arcTo(this.x, this.y + this._height, this.x, this.y + this._height - c, c), a.lineTo(this.x, this.y + c), a.arcTo(this.x, this.y, this.x + c, this.y, c), a.closePath(), null != this._alpha && (a.save(), a.globalAlpha = this._alpha), this.fill && (a.fillStyle = this.fill, a.fill()), this.stroke && (a.strokeStyle = this.stroke, a.lineWidth = 0.5, a.stroke())) : (null != this._alpha && (a.save(), a.globalAlpha = this._alpha), this.fill && (a.fillStyle = this.fill, a.fillRect(this.x, this.y, this._width, this._height)), this.stroke && (a.strokeStyle = this.stroke, a.lineWidth =
                0.5, a.strokeRect(this.x, this.y, this._width, this._height)));
            null != this._alpha && a.restore()
        };
        n.prototype.toSVG = function() {
            var a = x(A, "rect", null, {
                x: this.x,
                y: this.y,
                width: this._width,
                height: this._height,
                stroke: this.stroke || "none",
                strokeWidth: 0.5,
                fill: this.fill || "none"
            });
            null != this._alpha && a.setAttribute("opacity", this._alpha);
            return a
        };
        n.prototype.min = function() {
            return this.x
        };
        n.prototype.max = function() {
            return this.x + this._width
        };
        n.prototype.height = function() {
            return this.y + this._height
        };
        m.prototype.drawConnectors =
            function(a) {
                for (var c = this.coverage.ranges(), b = 1; b < c.length; ++b) {
                    var d = c[b],
                        f = c[b - 1];
                    if (f && d.min() > f.max()) {
                        var f = f.max(),
                            d = d.min(),
                            g = (f + d) / 2;
                        "hat+" === this.connector ? (a.moveTo(f, this.h / 2), a.lineTo(g, 0), a.lineTo(d, this.h / 2)) : "hat-" === this.connector ? (a.moveTo(f, this.h / 2), a.lineTo(g, this.h), a.lineTo(d, this.h / 2)) : "collapsed+" === this.connector ? (a.moveTo(f, this.h / 2), a.lineTo(d, this.h / 2), 4 < d - f && (a.moveTo(g - 2, this.h / 2 - 3), a.lineTo(g + 2, this.h / 2), a.lineTo(g - 2, this.h / 2 + 3))) : "collapsed-" === this.connector ? (a.moveTo(f,
                            this.h / 2), a.lineTo(d, this.h / 2), 4 < d - f && (a.moveTo(g + 2, this.h / 2 - 3), a.lineTo(g - 2, this.h / 2), a.lineTo(g + 2, this.h / 2 + 3))) : (a.moveTo(f, this.h / 2), a.lineTo(d, this.h / 2))
                    }
                }
            };
        m.prototype.draw = function(a, c) {
            for (var b = 0; b < this.glyphs.length; ++b) this.glyphs[b].draw(a, c);
            a.strokeStyle = "black";
            a.beginPath();
            this.drawConnectors(a);
            a.stroke()
        };
        m.prototype.toSVG = function() {
            for (var a = x(A, "g"), c = 0; c < this.glyphs.length; ++c) a.appendChild(this.glyphs[c].toSVG());
            c = new J;
            this.drawConnectors(c);
            0 < c.toPathData().length && (c = x(A,
                "path", null, {
                    d: c.toPathData(),
                    fill: "none",
                    stroke: "black",
                    strokeWidth: 0.5
                }), a.appendChild(c));
            return a
        };
        m.prototype.min = function() {
            return this.coverage.min()
        };
        m.prototype.max = function() {
            return this.coverage.max()
        };
        m.prototype.height = function() {
            return this.h
        };
        h.prototype.min = function() {
            return this.points[0]
        };
        h.prototype.max = function() {
            return this.points[this.points.length - 2]
        };
        h.prototype.height = function() {
            return this._height
        };
        h.prototype.draw = function(a) {
            a.save();
            a.strokeStyle = this.color;
            a.lineWidth =
                2;
            a.beginPath();
            for (var c = 0; c < this.points.length; c += 2) {
                var b = this.points[c],
                    d = this.points[c + 1];
                0 == c ? a.moveTo(b, d) : a.lineTo(b, d)
            }
            a.stroke();
            a.restore()
        };
        h.prototype.toSVG = function() {
            for (var a = new J, c = 0; c < this.points.length; c += 2) {
                var b = this.points[c],
                    d = this.points[c + 1];
                0 == c ? a.moveTo(b, d) : a.lineTo(b, d)
            }
            return x(A, "path", null, {
                d: a.toPathData(),
                fill: "none",
                stroke: this.color,
                strokeWidth: "2px"
            })
        };
        v.prototype.toSVG = function() {
            var a = this.glyph.toSVG(),
                c = {};
            "above" == this.align ? (a = x(A, "g", a, {
                transform: "translate(0, " +
                    (this.textHeight | 2) + ")"
            }), c.y = this.textHeight) : c.y = this.glyph.height() + 15;
            this.font && (c.fontSize = 7);
            c.x = "center" == this.anchor ? (this.glyph.min() + this.glyph.max() - this.textLen) / 2 : this.glyph.min();
            return x(A, "g", [a, x(A, "text", this.text, c)])
        };
        v.prototype.min = function() {
            return this.glyph.min()
        };
        v.prototype.max = function() {
            return this.measured ? Math.max(this.glyph.max(), 1 * this.glyph.min() + this.textLen + 10) : this.glyph.max()
        };
        v.prototype.height = function() {
            var a = this.glyph.height();
            this.measured && (a = "above" ==
                this.align ? a + (this.textHeight + 2) : a + 20);
            return a
        };
        v.prototype.draw = function(a, c) {
            "above" == this.align && (a.save(), a.translate(0, this.textHeight + 2));
            this.glyph.draw(a);
            "above" == this.align && a.restore();
            c.registerGlyph(this)
        };
        v.prototype.drawOverlay = function(a, c, b) {
            a.fillStyle = "black";
            this.font && (a.save(), a.font = this.font);
            "center" == this.anchor ? b = (this.glyph.min() + this.glyph.max() - this.textLen) / 2 : (b = this.glyph.min(), b < c && (b = Math.min(c, this.glyph.max() - this.textLen)));
            a.fillText(this.text, b, "above" == this.align ?
                this.textHeight : this.glyph.height() + 15);
            this.font && a.restore()
        };
        r.prototype.draw = function(a) {
            var c = this._height / 2;
            a.beginPath();
            a.moveTo(this._x, 0);
            a.lineTo(this._x, this._height);
            a.moveTo(this._x - c, c);
            a.lineTo(this._x + c, c);
            a.strokeStyle = this._stroke;
            a.lineWidth = 1;
            a.stroke()
        };
        r.prototype.toSVG = function() {
            var a = this._height / 2,
                c = new J;
            c.moveTo(this._x, 0);
            c.lineTo(this._x, this._height);
            c.moveTo(this._x - a, a);
            c.lineTo(this._x + a, a);
            return x(A, "path", null, {
                d: c.toPathData(),
                fill: "none",
                stroke: this._stroke,
                strokeWidth: "1px"
            })
        };
        r.prototype.min = function() {
            return this._x - this._height / 2
        };
        r.prototype.max = function() {
            return this._x + this._height / 2
        };
        r.prototype.height = function() {
            return this._height
        };
        z.prototype.draw = function(a) {
            var c = this._height / 2;
            a.beginPath();
            a.moveTo(this._x - c, 0);
            a.lineTo(this._x + c, this._height);
            a.moveTo(this._x - c, this._height);
            a.lineTo(this._x + c, 0);
            a.strokeStyle = this._stroke;
            a.lineWidth = 1;
            a.stroke()
        };
        z.prototype.toSVG = function() {
            var a = this._height / 2,
                c = new J;
            c.moveTo(this._x - a, 0);
            c.lineTo(this._x +
                a, this._height);
            c.moveTo(this._x - a, this._height);
            c.lineTo(this._x + a, 0);
            return x(A, "path", null, {
                d: c.toPathData(),
                fill: "none",
                stroke: this._stroke,
                strokeWidth: "1px"
            })
        };
        z.prototype.min = function() {
            return this._x - this._height / 2
        };
        z.prototype.max = function() {
            return this._x + this._height / 2
        };
        z.prototype.height = function() {
            return this._height
        };
        a.prototype = Object.create(p.prototype);
        a.prototype.drawPath = function(a) {
            var c = this._height / 2,
                b = this._width / 2;
            "S" === this._dir ? (a.moveTo(this._x, this._height), a.lineTo(this._x -
                b, 0), a.lineTo(this._x + b, 0)) : "W" === this._dir ? (a.moveTo(this._x + b, c), a.lineTo(this._x - b, 0), a.lineTo(this._x - b, this._height)) : "E" === this._dir ? (a.moveTo(this._x - b, c), a.lineTo(this._x + b, 0), a.lineTo(this._x + b, this._height)) : (a.moveTo(this._x, 0), a.lineTo(this._x + b, this._height), a.lineTo(this._x - b, this._height));
            a.closePath()
        };
        a.prototype.min = function() {
            return this._x - this._height / 2
        };
        a.prototype.max = function() {
            return this._x + this._height / 2
        };
        a.prototype.height = function() {
            return this._height
        };
        b.prototype.draw =
            function(a) {
                var c = this._height / 2;
                a.fillStyle = this._stroke;
                a.beginPath();
                a.arc(this._x, c, c, 0, 6.29);
                this._fill && (a.fillStyle = this._fill, a.fill());
                this._stroke && (a.strokeStyle = this._stroke, a.stroke())
            };
        b.prototype.toSVG = function() {
            var a = this._height / 2;
            return x(A, "circle", null, {
                cx: this._x,
                cy: a,
                r: a,
                fill: this._fill || "none",
                stroke: this._stroke || "none",
                strokeWidth: "1px"
            })
        };
        b.prototype.min = function() {
            return this._x - this._height / 2
        };
        b.prototype.max = function() {
            return this._x + this._height / 2
        };
        b.prototype.height =
            function() {
                return this._height
            };
        q.prototype.draw = function(a, c) {
            this.glyph && this.glyph.draw(a, c)
        };
        q.prototype.toSVG = function() {
            return this.glyph ? this.glyph.toSVG() : x(A, "g")
        };
        q.prototype.min = function() {
            return this._min
        };
        q.prototype.max = function() {
            return this._max
        };
        q.prototype.height = function() {
            return this.glyph ? this.glyph.height() : 1
        };
        g.prototype = Object.create(p.prototype);
        g.prototype.min = function() {
            return this._min
        };
        g.prototype.max = function() {
            return this._max
        };
        g.prototype.height = function() {
            return this._height
        };
        g.prototype.drawPath = function(a) {
            var c = this._max,
                b = this._min,
                d = this._height,
                f = 0,
                g = 0,
                l = this._height + 2,
                h = 0.333333 * this._height;
            this._ori && ("+" === this._ori ? g = 0.5 * this._height : "-" === this._ori && (f = 0.5 * this._height));
            c - b < l && (b = (c + b - l) / 2, c = b + l);
            a.moveTo(b + f, 0 + h);
            a.lineTo(c - g, 0 + h);
            a.lineTo(c - g, 0);
            a.lineTo(c, 0 + this._height / 2);
            a.lineTo(c - g, 0 + d);
            a.lineTo(c - g, 0 + h + h);
            a.lineTo(b + f, 0 + h + h);
            a.lineTo(b + f, 0 + d);
            a.lineTo(b, 0 + d / 2);
            a.lineTo(b + f, 0);
            a.lineTo(b + f, 0 + h)
        };
        l.prototype = Object.create(p.prototype);
        l.prototype.min =
            function() {
                return this._min
            };
        l.prototype.max = function() {
            return this._max
        };
        l.prototype.height = function() {
            return this._height
        };
        l.prototype.drawPath = function(a) {
            var c = this._min,
                b = this._max,
                d = this._height,
                f = d / 2;
            a.moveTo(c, f);
            a.lineTo(b, f);
            a.moveTo(c, 0);
            a.lineTo(c, d);
            a.moveTo(b, 0);
            a.lineTo(b, d)
        };
        f.prototype.min = function() {
            return this._min
        };
        f.prototype.max = function() {
            return this._max
        };
        f.prototype.height = function() {
            return this._height
        };
        f.prototype.drawPath = function(a) {
            var c = this._min,
                b = this._max,
                d = this._height,
                f = d / 2;
            "hat" === this._style ? (a.moveTo(c, f), a.lineTo((c + b) / 2, "-" === this._strand ? d : 0)) : a.moveTo(c, f);
            a.lineTo(b, f)
        };
        f.prototype.draw = function(a) {
            a.beginPath();
            this.drawPath(a);
            a.strokeStyle = this._stroke;
            "dashed" === this._style && a.setLineDash ? (a.save(), a.setLineDash([3]), a.stroke(), a.restore()) : a.stroke()
        };
        f.prototype.toSVG = function() {
            var a = new J;
            this.drawPath(a);
            a = {
                d: a.toPathData(),
                stroke: this._stroke || "none"
            };
            "dashed" === this._style && (a.strokeDasharray = "3");
            return x(A, "path", null, a)
        };
        d.prototype.min = function() {
            return this._min
        };
        d.prototype.max = function() {
            return this._max
        };
        d.prototype.height = function() {
            return this._height
        };
        d.prototype.drawStemPath = function(a) {
            var c = this._max,
                b = this._height / 2;
            a.moveTo(this._min, b);
            a.lineTo(c, b)
        };
        d.prototype.drawTrigsPath = function(a) {
            var c = this._min,
                b = this._max,
                d = this._height,
                f = d / 2;
            a.moveTo(c, 0);
            a.lineTo(c + d, f);
            a.lineTo(c, d);
            a.lineTo(c, 0);
            a.moveTo(b, 0);
            a.lineTo(b - d, f);
            a.lineTo(b, d);
            a.lineTo(b, 0)
        };
        d.prototype.draw = function(a) {
            a.beginPath();
            this.drawStemPath(a);
            a.strokeStyle = this._stroke;
            a.stroke();
            a.beginPath();
            this.drawTrigsPath(a);
            a.fillStyle = this._fill;
            a.fill()
        };
        d.prototype.toSVG = function() {
            var a = new J;
            this.drawStemPath(a);
            var c = new J;
            this.drawTrigsPath(c);
            return x(A, "g", [x(A, "path", null, {
                d: a.toPathData(),
                stroke: this._stroke || "none"
            }), x(A, "path", null, {
                d: c.toPathData(),
                fill: this._fill || "none"
            })])
        };
        t.prototype = Object.create(p.prototype);
        t.prototype.min = function() {
            return this._min
        };
        t.prototype.max = function() {
            return this._max
        };
        t.prototype.height = function() {
            return this._height
        };
        t.prototype.drawPath =
            function(a) {
                var c = this._min,
                    b = this._max,
                    d = this._height;
                if (this._parallel) {
                    var f = d / 2,
                        g = 0.4 * d;
                    this._sw ? (a.moveTo(c + f, d - g), a.lineTo(c + f, d), a.lineTo(c, f), a.lineTo(c + f, 0), a.lineTo(c + f, g)) : (a.moveTo(c, d - g), a.lineTo(c, g));
                    this._ne ? (a.lineTo(b - f, g), a.lineTo(b - f, 0), a.lineTo(b, f), a.lineTo(b - f, d), a.lineTo(b - f, d - g)) : (a.lineTo(b, g), a.lineTo(b, d - g))
                } else {
                    var f = (c + b) / 2,
                        g = 0.4 * (b - c),
                        l = d / 3;
                    this._ne ? (a.moveTo(c + g, l), a.lineTo(c, l), a.lineTo(f, 0), a.lineTo(b, l), a.lineTo(b - g, l)) : (a.moveTo(c + g, 0), a.lineTo(b - g, 0));
                    this._sw ?
                        (a.lineTo(b - g, d - l), a.lineTo(b, d - l), a.lineTo(f, d), a.lineTo(c, d - l), a.lineTo(c + g, d - l)) : (a.lineTo(b - g, d), a.lineTo(c + g, d))
                }
                a.closePath()
            };
        D.prototype.min = function() {
            return this._min
        };
        D.prototype.max = function() {
            return this._max
        };
        D.prototype.height = function() {
            return this._height
        };
        D.prototype.toSVG = function() {
            return x(A, "rect", null, {
                x: this._min,
                y: 0,
                width: this._max - this._min,
                height: this._height,
                stroke: this._stroke || "none",
                fill: this._fill || "none"
            })
        };
        D.prototype.draw = function(a) {
            this._fill && (a.fillStyle = this._fill,
                a.fillRect(this._min, 0, this._max - this._min, this._height));
            if (this._stroke) {
                a.strokeStyle = this._stroke;
                a.strokeRect(this._min, 0, this._max - this._min, this._height);
                a.beginPath();
                for (var c = 2; c < this._height; c += 3) a.moveTo(this._min, c), a.lineTo(this._max, c);
                a.stroke()
            }
        };
        P.prototype.min = function() {
            return this._min
        };
        P.prototype.max = function() {
            return Math.max(this._max, this._min + this._textLen)
        };
        P.prototype.height = function() {
            return this._height
        };
        P.prototype.draw = function(a) {
            a.fillStyle = this._fill;
            a.fillText(this._string,
                this._min, this._height - 4)
        };
        P.prototype.toSVG = function() {
            return x(A, "text", this._string, {
                x: this._min,
                y: this._height - 4
            })
        };
        N.prototype.min = function() {
            return this._min
        };
        N.prototype.max = function() {
            return this._max
        };
        N.prototype.height = function() {
            return this._height
        };
        N.prototype.draw = function(a) {
            var c = this._seq,
                b = this._fill;
            if (c) {
                var d = (this._max - this._min + 1) / c.length,
                    f = (3 - this._readframe) % 3,
                    g = (c.length - f) % 3,
                    f = "+" == this._orientation ? f : g;
                0 < f && (a.fillStyle = b, a.fillRect(this._min, 0, d * f, this._height));
                for (; f <
                    c.length; f += 3) {
                    g = c.substr(f, 3).toUpperCase();
                    "-" == this._orientation && (g = L(g));
                    var l = g in C ? C[g] : "?",
                        b = 3 == g.length ? G(l, f, this._fill) : this._fill;
                    a.fillStyle = b;
                    a.fillRect(this._min + f * d, 0, d * g.length, this._height);
                    8 <= d && 3 == g.length && (a.fillStyle = "white", a.fillText(l, this._min + (f + 1) * d, this._height))
                }
            }
        };
        N.prototype.toSVG = function() {
            var a = x(A, "g"),
                c = this._seq,
                b = this._fill;
            if (!c) return a;
            var d = (this._max - this._min + 1) / c.length,
                f = (3 - this._readframe) % 3,
                g = (c.length - f) % 3,
                f = "+" == this._orientation ? f : g;
            for (0 < f && a.appendChild(x(A,
                    "rect", null, {
                        x: this._min,
                        y: 0,
                        width: d * f,
                        height: this._height,
                        fill: b
                    })); f < c.length; f += 3) {
                g = c.substr(f, 3).toUpperCase();
                "-" == this._orientation && (g = L(g));
                var l = g in C ? C[g] : "?",
                    b = 3 == g.length ? G(l, f, this._fill) : this._fill;
                a.appendChild(x(A, "rect", null, {
                    x: this._min + f * d,
                    y: 0,
                    width: d * g.length,
                    height: this._height,
                    fill: b
                }));
                8 <= d && 3 == g.length && a.appendChild(x(A, "text", l, {
                    x: this._min + (f + 1) * d,
                    y: this._height,
                    fill: "white"
                }))
            }
            return a
        };
        (function(a) {
            function c(a, b, d, f, g, k, l, h, e, q) {
                this.baseColors = a;
                this._strandColor =
                    b;
                this._min = d;
                this._max = f;
                this._height = g;
                this._seq = k;
                this._ref = l;
                this._scheme = h;
                this._quals = e;
                this._fillbg = q
            }
            var b = 1 < window.devicePixelRatio,
                d = {},
                f = /^[ACGT-]$/;
            c.prototype.min = function() {
                return this._min
            };
            c.prototype.max = function() {
                return this._max
            };
            c.prototype.height = function() {
                return this._height
            };
            c.prototype.alphaForQual = function(a) {
                return 0.1 + 0.9 * Math.max(0, Math.min(1 * a / 30, 1))
            };
            c.prototype.draw = function(a) {
                var c = this._seq,
                    g = this._ref,
                    k = "mismatch-all" === this._scheme,
                    l = c ? c.length : this._max - this._min +
                    1,
                    h = (this._max - this._min + 1) / l;
                "mismatch" !== this._scheme && "mismatch-all" !== this._scheme || 8 <= h || (a.fillStyle = this._strandColor, a.fillRect(this._min, this._height / 4, this._max - this._min, this._height / 2));
                for (var e = 0; e < l; ++e) {
                    var q = c ? c.substr(e, 1).toUpperCase() : "N";
                    if (f.test(q) || 8 <= h) {
                        var w = this.baseColors[q];
                        if (this._quals) {
                            var t = this._quals.charCodeAt(e) - 33,
                                n = a.globalAlpha;
                            a.globalAlpha = this.alphaForQual(t)
                        }
                        w || (t = g ? g.substr(e, 1).toUpperCase() : "N", w = "N" == q || "N" == t ? "gray" : this._strandColor, k && (q = t));
                        a.fillStyle =
                            w;
                        t = f.test(q);
                        !this._fillbg && 8 <= h && t || a.fillRect(this._min + e * h, 0, h, this._height);
                        if (8 <= h && t) {
                            var t = w + "_" + q,
                                m = d[t];
                            if (!m) {
                                m = document.createElement("canvas");
                                b ? (m.width = 16, m.height = 20) : (m.width = 8, m.height = 10);
                                var p = m.getContext("2d");
                                b && p.scale(2, 2);
                                p.fillStyle = this._fillbg ? "black" : w;
                                w = p.measureText(q).width;
                                p.fillText(q, 0.5 * (8 - w), 8);
                                d[t] = m
                            }
                            b ? a.drawImage(m, this._min + e * h + 0.5 * (h - 8), 0, 8, 10) : a.drawImage(m, this._min + e * h + 0.5 * (h - 8), 0)
                        }
                        this._quals && (a.globalAlpha = n)
                    }
                }
            };
            c.prototype.toSVG = function() {
                for (var a =
                        this._seq, c = this._ref, b = "mismatch-all" === this._scheme, d = (this._max - this._min + 1) / this._seq.length, g = x(A, "g"), k = 0; k < a.length; ++k) {
                    var l = a ? a.substr(k, 1).toUpperCase() : "N",
                        h = this.baseColors[l];
                    if (!h) {
                        var e = c ? c.substr(k, 1).toUpperCase() : "N",
                            h = "N" == l || "N" == e ? "gray" : this._strandColor;
                        b && (l = e)
                    }
                    e = 1;
                    this._quals && (e = this._quals.charCodeAt(k) - 33, e = this.alphaForQual(e));
                    var q = f.test(l);
                    !this._fillbg && 8 <= d && q || g.appendChild(x(A, "rect", null, {
                        x: this._min + k * d,
                        y: 0,
                        width: d,
                        height: this._height,
                        fill: h,
                        fillOpacity: e
                    }));
                    8 <= d && q && g.appendChild(x(A, "text", l, {
                        x: this._min + (0.5 + k) * d,
                        y: 8,
                        textAnchor: "middle",
                        fill: this._fillbg ? "black" : h,
                        fillOpacity: e
                    }))
                }
                return g
            };
            a.SequenceGlyph = c
        })(this);
        Q.prototype.height = function() {
            return this._height ? this._height : this.glyph.height() + this._y
        };
        Q.prototype.min = function() {
            return this.glyph.min() + this._x
        };
        Q.prototype.max = function() {
            return this.glyph.max() + this._x
        };
        Q.prototype.minY = function() {
            return this._y
        };
        Q.prototype.maxY = function() {
            return this._y + this.glyph.height()
        };
        Q.prototype.draw =
            function(a, c) {
                a.save();
                a.translate(this._x, this._y);
                this.glyph.draw(a, c);
                a.restore()
            };
        Q.prototype.toSVG = function() {
            var a = this.glyph.toSVG();
            a.setAttribute("transform", "translate(" + this._x + "," + this._y + ")");
            return a
        };
        y.prototype.min = function() {
            return this._x - 2
        };
        y.prototype.max = function() {
            return this._x + 2
        };
        y.prototype.height = function() {
            return this._height
        };
        y.prototype.draw = function(a) {
            a.save();
            a.globalAlpha = 0.3;
            a.fillStyle = this._fill;
            a.beginPath();
            a.arc(this._x, this._y, 1.5, 0, 6.29);
            a.fill();
            a.restore()
        };
        y.prototype.toSVG = function() {
            return x(A, "circle", null, {
                cx: this._x,
                cy: this._y,
                r: 2,
                fill: this._fill,
                stroke: "none"
            })
        };
        H.prototype.notSelectable = !0;
        H.prototype.min = function() {
            return -1E5
        };
        H.prototype.max = function() {
            return 1E5
        };
        H.prototype.height = function() {
            return this._height
        };
        H.prototype.draw = function(a) {
            a.save();
            a.strokeStyle = "black";
            a.lineWidth = 0.1;
            a.beginPath();
            for (var c = 0; c <= this._height; c += 10) a.moveTo(-5E3, c), a.lineTo(5E3, c);
            a.stroke();
            a.restore()
        };
        H.prototype.toSVG = function() {
            for (var a = new J, c = 0; c <=
                this._height; c += 10) a.moveTo(-5E3, c), a.lineTo(5E3, c);
            return x(A, "path", null, {
                d: a.toPathData(),
                fill: "none",
                stroke: "black",
                strokeWidth: "0.1px"
            })
        };
        F.prototype = Object.create(p.prototype);
        F.prototype.min = function() {
            return this._x - this._r
        };
        F.prototype.max = function() {
            return this._x + this._r
        };
        F.prototype.height = function() {
            return 2 * this._r
        };
        F.prototype.drawPath = function(a) {
            for (var c = this._x, b = this._r, d = this._r, f = 0; f < this._points; ++f) {
                var g = 6.28 * f / this._points,
                    l = c + d * Math.sin(g),
                    g = b - d * Math.cos(g);
                0 == f ? a.moveTo(l,
                    g) : a.lineTo(l, g);
                g = 6.28 * (f + 0.5) / this._points;
                l = c + 0.4 * d * Math.sin(g);
                g = b - 0.4 * d * Math.cos(g);
                a.lineTo(l, g)
            }
            a.closePath()
        };
        c.prototype.draw = function(a) {
            var c = this._height / 2;
            a.fillStyle = this._stroke;
            a.beginPath();
            a.arc(this._x, c, c - this._overhang, 0, 6.29);
            a.moveTo(this._x, 0);
            a.lineTo(this._x, this._height);
            this._fill && (a.fillStyle = this._fill, a.fill());
            this._stroke && (a.strokeStyle = this._stroke, a.stroke())
        };
        c.prototype.toSVG = function() {
            var a = this._hh;
            return x(A, "g", [x(A, "circle", null, {
                    cx: this._x,
                    cy: a,
                    r: a - this._overhang
                }),
                x(A, "line", null, {
                    x1: this._x,
                    y1: 0,
                    x2: this._x,
                    y2: this._height
                })
            ], {
                fill: this._fill || "none",
                stroke: this._stroke || "none",
                strokeWidth: "1px"
            })
        };
        c.prototype.min = function() {
            return this._x - this._hh
        };
        c.prototype.max = function() {
            return this._x + this._hh
        };
        c.prototype.height = function() {
            return this._height
        };
        w.prototype.translate = function(a, c) {
            this.ox += a;
            this.oy += c
        };
        w.prototype.registerGlyph = function(a) {
            this.glyphs.push({
                x: this.ox,
                y: this.oy,
                glyph: a
            })
        };
        w.prototype.draw = function(a, c, b) {
            for (var d = 0; d < this.glyphs.length; ++d) {
                var f =
                    this.glyphs[d];
                a.save();
                a.translate(f.x, f.y);
                f.glyph.drawOverlay(a, c, b);
                a.restore()
            }
        };
        "undefined" !== typeof u && (u.exports = {
            BoxGlyph: n,
            GroupGlyph: m,
            LineGraphGlyph: h,
            LabelledGlyph: v,
            CrossGlyph: r,
            ExGlyph: z,
            TriangleGlyph: a,
            DotGlyph: b,
            PaddedGlyph: q,
            AArrowGlyph: g,
            SpanGlyph: l,
            LineGlyph: f,
            PrimersGlyph: d,
            ArrowGlyph: t,
            TooManyGlyph: D,
            TextGlyph: P,
            SequenceGlyph: this.SequenceGlyph,
            AminoAcidGlyph: N,
            TranslatedGlyph: Q,
            GridGlyph: H,
            StarGlyph: F,
            PointGlyph: y,
            PlimsollGlyph: c,
            OverlayLabelCanvas: w
        })
    }, {
        "./spans": 35,
        "./svg-utils": 38,
        "./utils": 48
    }],
    22: [function(e, u, s) {
        function p(e, h) {
            this.base = e;
            this.query = h
        }
        if ("undefined" !== typeof e) var n = e("./das").DASFeature;
        p.prototype.features = function(e, h, p) {
            url = this.base + "/features/" + e.name;
            h = [];
            this.query && h.push(this.query);
            e.isBounded && (h.push("start=" + e.start), h.push("end=" + e.end));
            0 < h.length && (url = url + "?" + h.join("&"));
            var r = new XMLHttpRequest;
            r.onreadystatechange = function() {
                if (4 == r.readyState)
                    if (300 <= r.status) p(null, "Error code " + r.status);
                    else {
                        var h = JSON.parse(r.response).features,
                            a = [];
                        for (fi = 0; fi < h.length; ++fi) {
                            var b = h[fi],
                                q = new n;
                            q.segment = e;
                            q.min = (b.start | 0) + 1;
                            q.max = b.end | 0;
                            b.name && (q.label = b.name);
                            q.type = b.type || "unknown";
                            a.push(q)
                        }
                        p(a)
                    }
            };
            r.open("GET", url, !0);
            r.responseType = "text";
            r.send("")
        };
        "undefined" !== typeof u && (u.exports = {
            JBrowseStore: p
        })
    }, {
        "./das": 10
    }],
    23: [function(e, u, s) {
        function p() {
            var a = this;
            this.reqs = [];
            this.awaitedFeatures = {};
            this.requestsIssued = new Q(function(b, d) {
                a.notifyRequestsIssued = b
            })
        }

        function n(a, b, d, c, f, g, l) {
            this.chr = a;
            this.min = b;
            this.max = d;
            this.coverage =
                l;
            this.scale = c;
            this.features = f || [];
            this.status = g
        }

        function m(a, b, d, c, g, l) {
            this.tierMap = a;
            this.chr = b;
            this.min = d;
            this.max = c;
            this.scale = g;
            this.seqSource = l || new f;
            this.viewCount = 0;
            this.featureCache = {};
            this.latestViews = {}
        }

        function h(a, b, d) {
            for (var c = [], f = {}, g = 0; g < a.length; ++g) {
                var l = a[g];
                l.min && l.max ? l.groups && 0 < l.groups.length ? r(f, l.groups[0].id, l) : l.min <= d && l.max >= b && c.push(l) : c.push(l)
            }
            for (var h in f) {
                a = f[h];
                for (var e = 1E11, q = -1E11, g = 0; g < a.length; ++g) l = a[g], e = Math.min(e, l.min), q = Math.max(q, l.max);
                if (e <=
                    d || q >= b)
                    for (g = 0; g < a.length; ++g) c.push(a[g])
            }
            return c
        }
        if ("undefined" !== typeof e) {
            s = e("./utils");
            var v = s.Awaited,
                r = s.pusho;
            s = e("./sourceadapters");
            var z = s.MappedFeatureSource,
                a = s.CachingFeatureSource,
                b = s.BWGFeatureSource,
                q = s.RemoteBWGFeatureSource,
                g = s.BAMFeatureSource,
                l = s.RemoteBAMFeatureSource,
                f = s.DummySequenceSource,
                d = s.DummyFeatureSource,
                t = e("./overlay").OverlayFeatureSource;
            s = e("./spans");
            var D = s.Range,
                P = s.intersection;
            s = e("./sample");
            var G = s.downsample,
                L = s.getBaseCoverage,
                N = e("./das").DASSequence,
                Q = e("es6-promise").Promise
        }
        p.prototype.addRequest = function(a) {
            this.reqs.push(a)
        };
        p.prototype.abortAll = function() {
            for (var a = 0; a < this.reqs.length; ++a) this.reqs[a].abort()
        };
        n.prototype.toString = function() {
            return this.chr + ":" + this.min + ".." + this.max + ";scale=" + this.scale
        };
        m.prototype.cancel = function() {
            this.cancelled = !0
        };
        m.prototype.bestCacheOverlapping = function(a, b, d) {
            return (a = this.featureCache[this.tierMap[0]]) ? a : null
        };
        m.prototype.viewFeatures = function(a, b, d, c) {
            if (c != c) throw "viewFeatures called with silly scale";
            if (a != this.chr) throw "Can't extend Known Space to a new chromosome";
            1 > b && (b = 1);
            this.min = b;
            this.max = d;
            this.scale = c;
            this.pool && this.pool.abortAll();
            this.pool = new p;
            this.awaitedSeq = new v;
            this.seqWasFetched = !1;
            this.viewCount++;
            this.startFetchesForTiers(this.tierMap);
            this.pool.notifyRequestsIssued()
        };
        m.prototype.invalidate = function(a) {
            this.pool && (this.featureCache[a] = null, this.startFetchesForTiers([a]))
        };
        m.prototype.startFetchesForTiers = function(a) {
            for (var b = this, d = this.awaitedSeq, c = !1, f, g = 0; g < a.length; ++g) try {
                this.startFetchesFor(a[g],
                    d) && (c = !0)
            } catch (l) {
                f = a[g], f.currentFeatures = [], f.currentSequence = null, f.draw(), f.updateHeight(), f.updateStatus(l), console.log("Error fetching tier source"), console.log(l), f = l
            }
            if (c && !this.seqWasFetched) {
                this.seqWasFetched = !0;
                var h = this.min,
                    e = this.max;
                if (this.cs && this.cs.start <= h && this.cs.end >= e) return a = this.cs.start == h && this.cs.end == e ? this.cs : new N(this.cs.name, h, e, this.cs.alphabet, this.cs.seq.substring(h - this.cs.start, e + 1 - this.cs.start)), d.provide(a);
                this.seqSource.fetch(this.chr, h, e, this.pool,
                    function(a, c) {
                        if (c) {
                            if (!b.cs || h <= b.cs.start && e >= b.cs.end || h >= b.cs.end || e <= b.cs.start || e - h > b.cs.end - b.cs.start) b.cs = c;
                            d.provide(c)
                        } else console.log("Sequence loading failed", a), d.provide(null)
                    })
            }
            if (f) throw f;
        };
        m.prototype.startFetchesFor = function(a, b) {
            var f = this,
                c = this.viewCount,
                g = a.getSource() || new d,
                l = a.needsSequence(this.scale),
                e = f.featureCache[a],
                q = a.getActiveStyleFilters(this.scale),
                t;
            q && (t = q.typeList());
            var m = this.chr,
                p = this.min,
                k = this.max;
            if (void 0 === t) return !1;
            if (e && e.chr === this.chr && e.min <=
                p && e.max >= k) {
                var r = e.features;
                if (e.min < p || e.max > k) r = h(r, p, k);
                f.provision(a, e.chr, P(e.coverage, new D(p, k)), e.scale, t, r, e.status, l ? b : null);
                r = g.getScales();
                if (e.scale <= this.scale || !r) return l
            }
            g.instrument && console.log("Starting  fetch " + c + " (" + p + ", " + k + ")");
            g.fetch(m, p, k, this.scale, t, this.pool, function(d, h, q, r) {
                g.instrument && console.log("Finishing fetch " + c);
                var x = f.latestViews[a] || -1;
                if (!(f.cancelled || x > c)) {
                    r || (r = new D(p, k));
                    if (!e || p < e.min || k > e.max) f.featureCache[a] = new n(m, p, k, q, h, d, r);
                    f.latestViews[a] =
                        c;
                    f.provision(a, m, r, q, t, h, d, l ? b : null)
                }
            }, q);
            return l
        };
        m.prototype.provision = function(d, f, h, c, e, m, n, p) {
            n && (d.currentFeatures = [], d.currentSequence = null, d.draw(), d.updateHeight());
            d.updateStatus(n);
            if (!n) {
                for (var r = n = !1, A = d.getSource(); z.prototype.isPrototypeOf(A) || a.prototype.isPrototypeOf(A) || t.prototype.isPrototypeOf(A);) A = t.prototype.isPrototypeOf(A) ? A.sources[0] : A.source;
                if (b.prototype.isPrototypeOf(A) || q.prototype.isPrototypeOf(A) || g.prototype.isPrototypeOf(A) || l.prototype.isPrototypeOf(A)) n = !0;
                A.opts && (A.opts.forceReduction || A.opts.noDownsample) || n && e && 1 == e.length && 0 <= e.indexOf("density") && (m = G(m, this.scale));
                e && 1 == e.length && 0 <= e.indexOf("base-coverage") && (r = !0);
                p ? p.await(function(a) {
                    r && (m = L(m, a, d.browser.baseColors));
                    d.viewFeatures(f, h, c, m, a)
                }) : d.viewFeatures(f, h, c, m)
            }
        };
        "undefined" !== typeof u && (u.exports = {
            KnownSpace: m
        })
    }, {
        "./das": 10,
        "./overlay": 27,
        "./sample": 29,
        "./sourceadapters": 34,
        "./spans": 35,
        "./utils": 48,
        "es6-promise": 52
    }],
    24: [function(e, u, s) {
        function p(a, h) {
            this.block = a;
            this.offset =
                h
        }

        function n(a, h) {
            var g = 4294967296 * (a[h + 6] & 255) + 16777216 * (a[h + 5] & 255) + 65536 * (a[h + 4] & 255) + 256 * (a[h + 3] & 255) + (a[h + 2] & 255),
                l = a[h + 1] << 8 | a[h];
            return 0 == g && 0 == l ? null : new p(g, l)
        }

        function m(b, h) {
            h = Math.min(h || 1, b.byteLength - 50);
            for (var g = [], l = [0], f = 0; l[0] < h;) {
                var d = new Uint8Array(b, l[0], 12),
                    d = d[11] << 8 | d[10],
                    d = z(b, 12 + d + l[0], Math.min(65536, b.byteLength - 12 - d - l[0]), l);
                l[0] += 8;
                f += d.byteLength;
                g.push(d)
            }
            if (1 == g.length) return g[0];
            l = new Uint8Array(f);
            for (d = f = 0; d < g.length; ++d) {
                var e = new Uint8Array(g[d]);
                a(e, 0, l, f,
                    e.length);
                f += e.length
            }
            return l.buffer
        }

        function h(a, h) {
            this.minv = a;
            this.maxv = h
        }

        function v(a, h) {
            --h;
            return a >> 14 == h >> 14 ? 4681 + (a >> 14) : a >> 17 == h >> 17 ? 585 + (a >> 17) : a >> 20 == h >> 20 ? 73 + (a >> 20) : a >> 23 == h >> 23 ? 9 + (a >> 23) : a >> 26 == h >> 26 ? 1 + (a >> 26) : 0
        }

        function r(a, h) {
            var g, l = [];
            --h;
            l.push(0);
            for (g = 1 + (a >> 26); g <= 1 + (h >> 26); ++g) l.push(g);
            for (g = 9 + (a >> 23); g <= 9 + (h >> 23); ++g) l.push(g);
            for (g = 73 + (a >> 20); g <= 73 + (h >> 20); ++g) l.push(g);
            for (g = 585 + (a >> 17); g <= 585 + (h >> 17); ++g) l.push(g);
            for (g = 4681 + (a >> 14); g <= 4681 + (h >> 14); ++g) l.push(g);
            return l
        }
        if ("undefined" !== typeof e) {
            e = e("jszlib");
            var z = e.inflateBuffer,
                a = e.arrayCopy
        }
        p.prototype.toString = function() {
            return "" + this.block + ":" + this.offset
        };
        "undefined" !== typeof u && (u.exports = {
            unbgzf: m,
            readVob: n,
            reg2bin: v,
            reg2bins: r,
            Chunk: h
        })
    }, {
        jszlib: 64
    }],
    25: [function(e, u, s) {
        function p() {
            this.featuresByChr = {};
            this.maxLength = 1;
            this.chrRing = null
        }

        function n(a) {
            this.source = a;
            v.call(this);
            this.storeHolder = new r;
            this.parser = h(a.payload);
            if (!this.parser) throw "Unsupported memstore payload: " + a.payload;
            var b = this;
            this._load(function(a,
                g) {
                if (a) {
                    for (var l = new p, f = [], d = a.split("\n"), h = b.parser.createSession(function(a) {
                            f.push(a)
                        }), e = 0; e < d.length; ++e) {
                        var m = d[e];
                        0 < m.length && h.parse(m)
                    }
                    h.flush();
                    l.addFeatures(f);
                    b.storeHolder.provide(l)
                } else b.error = g || "No data", b.storeHolder.provide(null)
            })
        }
        if ("undefined" !== typeof e) {
            u = e("./sourceadapters");
            var m = u.registerSourceAdapterFactory,
                h = u.makeParser,
                v = u.FeatureSourceBase;
            e("./das");
            e = e("./utils");
            var r = e.Awaited,
                z = e.textXHR
        }
        p.prototype.addFeatures = function(a) {
            for (var b = {}, h = 0; h < a.length; ++h) {
                var g =
                    a[h],
                    l = g.segment || g.chr,
                    f = this.featuresByChr[l];
                f || (f = [], this.featuresByChr[l] = f);
                f.push(g);
                b[l] = !0;
                g = g.max - g.min + 1;
                g > this.maxLength && (this.maxLength = g)
            }
            for (l in b) f = this.featuresByChr[l], f.sort(function(a, b) {
                var f = a.min - b.min;
                return 0 != f ? f : a.max - b.max
            });
            this.chrRing = null
        };
        p.prototype._indexFor = function(a, b) {
            for (var h = 0, g = a.length; g > h;) {
                var l = (h + g) / 2 | 0;
                if (l >= a.length) return a.length;
                b < a[l].min ? g = l : h = l + 1
            }
            return g
        };
        p.prototype.fetch = function(a, b, h) {
            var g = this.featuresByChr[a];
            g || (g = 0 == a.indexOf("chr") ?
                this.featuresByChr[a.substring(3)] : this.featuresByChr["chr" + a]);
            if (!g) return [];
            var l = Math.max(0, this._indexFor(g, b - this.maxLength - 1));
            a = Math.min(g.length - 1, this._indexFor(g, h));
            for (var f = []; l <= a; ++l) {
                var d = g[l];
                d.min <= h && d.max >= b && f.push(d)
            }
            return f
        };
        p.prototype.findNextFeature = function(a, b, h) {
            if (null == this.chrRing) {
                this.chrRing = [];
                for (a in this.featuresByChr) this.chrRing.push(a);
                this.chrRing.sort()
            }
            var g = this.featuresByChr[a];
            g || (a = 0 == a.indexOf("chr") ? a.substring(3) : "chr" + a, g = this.featuresByChr[a]);
            if (!g) return null;
            var l = Math.max(0, Math.min(this._indexFor(g, b), g.length - 1));
            if (0 < h) {
                for (; l < g.length;) {
                    var f = g[l++];
                    if (f.min > b) return f
                }
                a = this.chrRing.indexOf(a) + 1;
                a >= this.chrRing.length && (a = 0);
                return this.findNextFeature(this.chrRing[a], 0, h)
            }
            for (; 0 <= l;)
                if (f = g[l--], f.max < b) return f;
            a = this.chrRing.indexOf(a) - 1;
            0 > a && (a = this.chrRing.length - 1);
            return this.findNextFeature(this.chrRing[a], 1E10, h)
        };
        n.prototype = Object.create(v.prototype);
        n.prototype._load = function(a) {
            if (this.source.blob) {
                var b = new FileReader;
                b.onloadend = function() {
                    return a(b.result, b.error)
                };
                b.readAsText(this.source.blob)
            } else {
                if (this.source.credentials) var h = {
                    credentials: this.source.credentials
                };
                z(this.source.uri, a, h)
            }
        };
        n.prototype.fetch = function(a, b, h, g, l, f, d) {
            var e = this;
            this.storeHolder.await(function(f) {
                return f ? (f = f.fetch(a, b, h), d(null, f, 1E8)) : d(e.error)
            })
        };
        n.prototype.getStyleSheet = function(a) {
            this.parser && this.parser.getStyleSheet && this.parser.getStyleSheet(a)
        };
        n.prototype.getDefaultFIPs = function(a) {
            this.parser && this.parser.getDefaultFIPs &&
                this.parser.getDefaultFIPs(a)
        };
        n.prototype.getScales = function() {
            return 1E8
        };
        n.prototype.findNextFeature = function(a, b, h, g) {
            var l = this;
            this.storeHolder.await(function(f) {
                return f ? g(f.findNextFeature(a, b, h)) : g(null, l.error)
            })
        };
        n.prototype.capabilities = function() {
            return {
                leap: !0
            }
        };
        m("memstore", function(a) {
            return {
                features: new n(a)
            }
        })
    }, {
        "./das": 10,
        "./sourceadapters": 34,
        "./utils": 48
    }],
    26: [function(e, u, s) {
        function p(e) {
            return (e | 0).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",")
        }

        function n(e) {
            e = "" + e;
            var h =
                e.indexOf(".");
            if (0 > h) return e;
            var n = 2;
            "-" == e.substring(0, 1) && ++n;
            return h >= n ? e.substring(0, h) : e.substring(0, h + 2)
        }
        "undefined" !== typeof u && (u.exports = {
            formatLongInt: p,
            formatQuantLabel: n
        })
    }, {}],
    27: [function(e, u, s) {
        function p(e, n) {
            this.sources = e;
            this.opts = n || {};
            this.activityListeners = [];
            this.readinessListeners = [];
            this.changeListeners = [];
            this.business = [];
            this.readiness = [];
            for (var a = 0; a < this.sources.length; ++a) this.initN(a);
            "function" === typeof n.merge ? this.merge = n.merge : "concat" == n.merge ? this.merge = h :
                "alternates" == n.merge ? (this.merge = h, this.filterDispatchOnMethod = !0) : this.merge = m
        }

        function n(h, e, a) {
            this.source = h;
            this.callback = e;
            this.sources = a;
            this.count = a.length;
            this.statusCount = this.returnCount = 0;
            this.returns = [];
            this.features = [];
            this.statuses = [];
            this.scale = null
        }

        function m(h) {
            for (var e = [], a = 1; a < h.length; ++a) {
                for (var b = {}, q = h[a], g = 0; g < q.length; ++g) b[this.keyForFeature(q[g])] = q[g];
                e.push(b)
            }
            a = [];
            h = h[0];
            for (g = 0; g < h.length; ++g) {
                for (var l = h[g], f = 0; f < e.length; ++f)
                    if (b = e[f], q = b[this.keyForFeature(l)])
                        for (var d in q) "score" ===
                            d ? l.score2 = q.score : "min" !== d && "max" !== d && "segment" !== d && "_cachedStyle" !== d && (l[d] = q[d]);
                a.push(l)
            }
            return a
        }

        function h(h, e) {
            for (var a = [], b = 0; b < h.length; ++b)
                for (var q = h[b], g = e[b].name, l = 0; l < q.length; ++l) {
                    var f = q[l];
                    f.method = g;
                    a.push(f)
                }
            return a
        }
        if ("undefined" !== typeof e) var v = e("./utils").shallowCopy;
        p.prototype.initN = function(h) {
            var e = this.sources[h],
                a = this;
            this.business[h] = 0;
            e.addActivityListener && e.addActivityListener(function(b) {
                a.business[h] = b;
                a.notifyActivity()
            });
            e.addChangeListener && e.addChangeListener(function() {
                a.notifyChange()
            });
            e.addReadinessListener && e.addReadinessListener(function(b) {
                a.readiness[h] = b;
                a.notifyReadiness()
            })
        };
        p.prototype.addReadinessListener = function(h) {
            this.readinessListeners.push(h);
            this.notifyReadinessListener(h)
        };
        p.prototype.notifyReadiness = function() {
            for (var h = 0; h < this.readinessListeners.length; ++h) this.notifyReadinessListener(this.readinessListeners[h])
        };
        p.prototype.notifyReadinessListener = function(h) {
            for (var e = null, a = 0; a < this.readiness.length; ++a)
                if (null != this.readiness[a]) {
                    e = this.readiness[a];
                    break
                }
            try {
                h(e)
            } catch (b) {
                console.log(b)
            }
        };
        p.prototype.addActivityListener = function(h) {
            this.activityListeners.push(h)
        };
        p.prototype.notifyActivity = function() {
            for (var h = 0, e = 0; e < this.business.length; ++e) h += this.business[e];
            for (e = 0; e < this.activityListeners.length; ++e) try {
                this.activityListeners[e](h)
            } catch (a) {
                console.log(a)
            }
        };
        p.prototype.addChangeListener = function(h) {
            this.changeListeners.push(h)
        };
        p.prototype.notifyChange = function() {
            for (var h = 0; h < this.changeListeners.length; ++h) try {
                this.changeListeners[h](this.busy)
            } catch (e) {
                console.log(e)
            }
        };
        p.prototype.getScales =
            function() {
                return this.sources[0].getScales()
            };
        p.prototype.getStyleSheet = function(h) {
            return this.sources[0].getStyleSheet(h)
        };
        p.prototype.capabilities = function() {
            var h = {},
                e = this.sources[0];
            e.capabilities && (h = v(e.capabilities()));
            for (e = 1; e < this.sources.length; ++e) {
                var a = this.sources[e];
                a.capabilities && (a = a.capabilities(), a.search && (h.search = a.search))
            }
            return h
        };
        p.prototype.search = function(h, e) {
            for (var a = 0; a < this.sources.length; ++a) {
                var b;
                b = this.sources[a];
                b = b.capabilities ? b.capabilities().search : !1;
                if (b) return this.sources[a].search(h, e)
            }
        };
        p.prototype.fetch = function(h, e, a, b, q, g, l, f) {
            var d;
            if (this.filterDispatchOnMethod) {
                d = [];
                for (var t = f.list(), m = 0; m < this.sources.length; ++m)
                    for (var p = this.sources[m], G = 0; G < t.length; ++G) {
                        var L = t[G];
                        if (!L.method || L.method == p.name) {
                            d.push(p);
                            break
                        }
                    }
            } else d = this.sources;
            l = new n(this, l, d);
            for (m = 0; m < d.length; ++m) this.fetchN(l, m, d[m], h, e, a, b, q, g, f)
        };
        p.prototype.fetchN = function(h, e, a, b, q, g, l, f, d, n) {
            a.fetch(b, q, g, l, f, d, function(a, b, d) {
                return h.completed(e, a, b, d)
            }, n)
        };
        p.prototype.quantFindNextFeature = function(h, e, a, b, q) {
            return this.sources[0].quantFindNextFeature(h, e, a, b, q)
        };
        p.prototype.findNextFeature = function(h, e, a, b) {
            return this.sources[0].findNextFeature(h, e, a, b)
        };
        n.prototype.completed = function(h, e, a, b) {
            if (null == this.scale || 0 == h) this.scale = b;
            if (this.returns[h]) throw "Multiple returns for source " + h;
            this.returns[h] = !0;
            this.returnCount++;
            this.features[h] = a;
            e && (this.statuses[h] = e, this.statusCount++);
            if (this.returnCount == this.count) {
                if (0 < this.statusCount) {
                    h = "";
                    for (e = 0; e < this.count; ++e)
                        if (a = this.statuses[e]) 0 < h.length && (h += ", "), h += a;
                    return this.callback(h, null, this.scale)
                }
                this.callback(null, this.source.merge(this.features, this.sources), this.scale)
            }
        };
        p.prototype.getDefaultFIPs = function(h) {
            for (var e = 0; e < this.sources.length; ++e) {
                var a = this.sources[e];
                a.getDefaultFIPs && a.getDefaultFIPs(h)
            }
        };
        p.prototype.keyForFeature = function(h) {
            return "" + h.min + ".." + h.max
        };
        "undefined" !== typeof u && (u.exports = {
            OverlayFeatureSource: p
        })
    }, {
        "./utils": 48
    }],
    28: [function(e, u, s) {
        function p(f,
            d, e) {
            var D = /^\w+\s[0-9]+\s[0-9]+.*$/,
                s = /([^=]+)=\"?([^\"]+)\"?/,
                G = /^##\s*fileformat=VCFv4\..+/;
            (f.blob ? new h(f.blob) : "encode" == f.transport ? new l(f.uri) : new m(f.uri, {
                credentials: f.credentials
            })).slice(0, 65536).salted().fetch(function(l, h) {
                if (!l) return e || (f.credentials = !0, p(f, d, !0)), d(f, "Couldn't fetch data");
                var m = new Uint8Array(l),
                    y = (new Uint32Array(l, 0, 1))[0];
                if (y == r || y == z) {
                    f.tier_type = "bwg";
                    m = /\/?([^/]+?)(.bw|.bb|.bigWig|.bigBed)?$/;
                    if (m = m.exec(f.uri || f.blob.name)) f.name = m[1];
                    return d(f, null)
                }
                if (y ==
                    q) return f.tier_type = "bai", d(f, null);
                if (31 == m[0] || 139 == m[1]) {
                    m = a(l);
                    m = new Uint8Array(m);
                    y = v(m, 0);
                    if (y == b) {
                        f.tier_type = "bam";
                        m = /\/?([^/]+?)(.bam)?$/;
                        if (m = m.exec(f.uri || f.blob.name)) f.name = m[1];
                        return d(f, null)
                    }
                    if (y == g) return f.tier_type = "tabix-index", d(f, null);
                    if (1768301347 == y) {
                        f.tier_type = "tabix";
                        f.payload = "vcf";
                        m = /\/?([^/]+?)(.vcf)?(.gz)?$/;
                        if (m = m.exec(f.uri || f.blob.name)) f.name = m[1];
                        return d(f, null)
                    }
                    console.log("magic = " + y.toString(16));
                    return d(f, "Unsupported format")
                }
                m = String.fromCharCode.apply(null,
                    m).split("\n");
                if (0 < m.length && G.test(m[0])) return f.tier_type = "memstore", f.payload = "vcf", m = /\/?([^/]+?)(.vcf)?$/, (m = m.exec(f.uri || f.blob.name)) && !f.name && (f.name = m[1]), d(f, null);
                for (var H = 0; H < m.length; ++H)
                    if (y = m[H].replace("\r", ""), 0 != y.length && 0 != y.indexOf("browser")) {
                        if (0 == y.indexOf("track")) {
                            m = "bed";
                            y = y.split(/\s/);
                            for (H = 1; H < y.length; ++H) {
                                var u = s.exec(y[H]);
                                u && ("type" == u[1] && "wiggle_0" == u[2] ? m = "wig" : "name" == u[0] && (f.name = u[2]))
                            }
                            n(f, m);
                            return d(f, null)
                        }
                        if (0 == y.indexOf("fixedStep") || 0 == y.indexOf("variableStep")) return n(f,
                            "wig"), d(f, null);
                        if (D.test(y)) return n(f, null), d(f, null);
                        break
                    }
                return d(f, "Unsupported format")
            })
        }

        function n(a, b) {
            a.tier_type = "memstore";
            var g = /\/?([^/]+?)(.(bed|wig))?$/.exec(a.uri || a.blob.name);
            g && (a.name || (a.name = g[1]), !b && g[3] && (b = g[3]));
            a.payload = b || "bed"
        }
        if ("undefined" !== typeof e) {
            s = e("./bin");
            var m = s.URLFetchable,
                h = s.BlobFetchable,
                v = s.readInt;
            s = e("./bigwig");
            var r = s.BIG_WIG_MAGIC,
                z = s.BIG_BED_MAGIC,
                a = e("./lh3utils").unbgzf;
            s = e("./bam");
            var b = s.BAM_MAGIC,
                q = s.BAI_MAGIC,
                g = e("./tabix").TABIX_MAGIC,
                l = e("./encode").EncodeFetchable
        }
        "undefined" !== typeof u && (u.exports = {
            probeResource: p
        })
    }, {
        "./bam": 1,
        "./bigwig": 3,
        "./bin": 4,
        "./encode": 12,
        "./lh3utils": 24,
        "./tabix": 40
    }],
    29: [function(e, u, s) {
        function p(a) {
            return q[a % q.length] * Math.pow(10, a / q.length | 0)
        }

        function n(a, b, f) {
            this.scale = a;
            this.cnt = this.tot = 0;
            this.hasScore = !1;
            this.min = b;
            this.max = f;
            this.features = []
        }

        function m(a, b) {
            return a.min < b.min ? -1 : a.min > b.min ? 1 : a.max < b.max ? -1 : b.max > a.max ? 1 : 0
        }

        function h(a, b) {
            for (var f = 0; p(f + 1) < b;) ++f;
            for (var f = p(f), d = [], h = -1E10, e = 1E10, q = 0; q < a.length; ++q) {
                var m = a[q];
                if (m.groups && 0 < m.groups.length) return a;
                for (var r = m.min / f | 0, s = m.max / f | 0, h = Math.max(h, s), e = Math.min(e, r); r <= s; ++r) {
                    var v = d[r];
                    v || (v = new n(f, r * f, (r + 1) * f - 1), d[r] = v);
                    v.feature(m)
                }
            }
            q = [];
            for (r = e; r <= h; ++r)
                if (v = d[r]) m = new z, m.segment = a[0].segment, m.min = r * f + 1, m.max = (r + 1) * f, m.score = v.score(), m.type = "density", q.push(m);
            Date.now();
            return q
        }

        function v(a) {
            this._pos = a;
            this._bases = {};
            this._totalCount = 0
        }

        function r(g, l, f) {
            for (var d, h, e = null, q = null, m = [], n = 0; n < g.length; ++n) {
                var p =
                    g[n];
                if (p.groups && 0 < p.groups.length) return g;
                d = p.seq;
                var r = p.quals,
                    y = a(p.cigar),
                    s = [];
                h = [];
                for (var u = 0, c = 0; c < y.length; ++c) {
                    var w = y[c];
                    if ("M" == w.op) s.push(d.substr(u, w.cnt)), h.push(r.substr(u, w.cnt)), u += w.cnt;
                    else if ("D" == w.op)
                        for (var B = 0; B < w.cnt; ++B) s.push("-"), h.push("Z");
                    else "I" == w.op ? u += w.cnt : "S" == w.op ? u += w.cnt : console.log("unknown cigop" + w.op)
                }
                d = s.join("");
                h = h.join("");
                s = d;
                u = h;
                c = p.orientation;
                w = p.min || 0;
                p = p.max || 0;
                h = 0;
                for (d = w; d <= p; ++d) y = m[d], y || (y = new v(d), m[d] = y), r = s.charAt(h), B = u.charCodeAt(h) -
                    33, y.recordBase(r, B, c), h++;
                e = e ? Math.min(e, w) : w;
                q = q ? Math.max(q, p) : p
            }
            p = e;
            n = q;
            h = [];
            if (l && (d = l.start | 0, y = l.end | 0, d <= n && y >= p)) {
                r = Math.max(p, d);
                y = Math.min(n, y);
                for (s = 0; s < r - p; s++) h.push("N");
                h.push(l.seq.substr(r - d, y - r + 1));
                for (s = 0; s < n - y; s++) h.push("N")
            }
            l = h.join("");
            n = [];
            h = 0;
            for (d = e; d <= q; ++d) {
                if (y = m[d])
                    if (p = new z, p.segment = g[0].segment, p.min = y.pos(), p.max = p.min, p.notes = [], p.notes = p.notes.concat(y.infoList()), p.type = "base-coverage", p.suppressScore = !0, l)
                        for (e = l.charAt(h), p.notes.unshift("Ref=" + e), y = y.baseScoreList(e,
                                0.2), s = 0; s < y.length; s++) {
                            r = y[s].base;
                            u = y[s].score;
                            c = b(p);
                            c.score = u;
                            if (1 < y.length || r != e) c.itemRgb = f[r];
                            n.push(c)
                        } else n.push(p);
                h++
            }
            return n
        }
        if ("undefined" !== typeof e) var z = e("./das").DASFeature,
            a = e("./cigar").parseCigar,
            b = e("./utils").shallowCopy;
        var q = [1, 2, 5];
        n.prototype.score = function() {
            if (0 == this.cnt) return 0;
            if (this.hasScore) return this.tot / this.cnt;
            var a = this.features;
            a.sort(m);
            for (var b = -1E10, f = 0, d = 0, h = 1; h < a.length; ++h) {
                var e = a[h],
                    q = Math.max(e.min, this.min),
                    e = Math.min(e.max, this.max),
                    d = d + (e -
                        q + 1);
                q > b ? (f += e - q + 1, b = e) : e > b && (f += e - b, b = e)
            }
            return 0 < f ? 1 * d / f : 0
        };
        n.prototype.feature = function(a) {
            void 0 !== a.score && (this.tot += a.score, this.hasScore = !0);
            ++this.cnt;
            this.features.push(a)
        };
        v.prototype.recordBase = function(a, b, f) {
            if (this._bases[a]) a = this._bases[a], a.cnt++, a.totalQual += b, a.strandCnt[f]++;
            else {
                var d = {
                    "+": 0,
                    "-": 0
                };
                d[f]++;
                this._bases[a] = {
                    cnt: 1,
                    totalQual: b,
                    strandCnt: d
                }
            }
            this._totalCount++
        };
        v.prototype.totalCount = function() {
            return this._totalCount
        };
        v.prototype.pos = function() {
            return this._pos
        };
        v.prototype.infoList =
            function() {
                var a = [],
                    b = this._totalCount;
                a.push("Depth=" + b.toString());
                for (var f in this._bases) {
                    var d = this._bases[f],
                        h = d.cnt,
                        e = d.strandCnt["+"],
                        q = d.strandCnt["-"],
                        d = d.totalQual / h,
                        h = [f, "=", h, " (", (100 * h / b).toFixed(0), "%, ", e, " +, ", q, " -, Qual: ", d.toFixed(0), ")"];
                    a.push(h.join(""))
                }
                return a
            };
        v.prototype.baseScoreList = function(a, b) {
            var f = [],
                d = this._totalCount,
                h = b * d,
                e;
            for (e in this._bases) {
                var q = this._bases[e].cnt;
                q < h || e == a || (f.push({
                    base: e,
                    score: d
                }), d -= q)
            }
            f.push({
                base: a,
                score: d
            });
            return f
        };
        "undefined" !==
        typeof u && (u.exports = {
            downsample: h,
            getBaseCoverage: r
        })
    }, {
        "./cigar": 8,
        "./das": 10,
        "./utils": 48
    }],
    30: [function(e, u, s) {
        function p(h, e) {
            var a = parseFloat(h.replace(/,/g, ""));
            return "k" === e || "K" === e ? 1E3 * a | 0 : "m" == e || "M" === e ? 1E6 * a | 0 : a | 0
        }
        if ("undefined" !== typeof e) var n = e("./cbrowser").Browser,
            m = e("./bin").URLFetchable,
            h = e("./trix").connectTrix;
        var v = /^([\d+,\w,\.,\_,\-]+):([0-9,\.]+?)([KkMmGg])?((-|\.\.)+([0-9,\.]+)([KkMmGg])?)?$/;
        n.prototype.search = function(e, n) {
            var a = this,
                b = v.exec(e);
            if (b) {
                var q = b[1],
                    g;
                if (b[6]) g =
                    p(b[2], b[3]), b = p(b[6], b[7]);
                else {
                    var l = this.viewEnd - this.viewStart + 1;
                    g = p(b[2], b[3]) - l / 2 | 0;
                    b = g + l - 1
                }
                this.setLocation(q, g, b, n)
            } else {
                if (!e || 0 == e.length) return !1;
                var f = 0,
                    d = !1,
                    t = function(b, g) {
                        --f;
                        if (g) return n(g);
                        b || (b = []);
                        for (var h = 5E8, l = -1E8, q = null, m = 0; m < b.length; ++m) {
                            var p = b[m];
                            null == q && (q = p.segment);
                            h = Math.min(h, p.min);
                            l = Math.max(l, p.max)
                        }
                        if (q) d = !0, a.highlightRegion(q, h, l), m = Math.max(2500, 0.3 * (l - h + 1) | 0), a.setLocation(q, h - m, l + m, n);
                        else if (0 == f && !d) return n("no match for '" + e + "'")
                    },
                    D = function(a, b) {
                        b.lookup(e,
                            function(b, d) {
                                if (null == b || 2 > b.length) return a.featureSource.search(e, t);
                                var f = b[1].split(",")[0];
                                return a.featureSource.search(f, t)
                            })
                    };
                if (this.searchEndpoint) return f = 1, this.doDasSearch(a.searchEndpoint, e, t);
                for (q = 0; q < this.tiers.length; ++q)(function(b) {
                    a.sourceAdapterIsCapable(b.featureSource, "search") ? b.dasSource.trixURI ? (++f, b.trix ? D(b, b.trix) : h(new m(b.dasSource.trixURI), new m(b.dasSource.trixURI + "x"), function(a) {
                            b.trix = a;
                            D(b, a)
                        })) : (++f, b.featureSource.search(e, t)) : b.dasSource.provides_search &&
                        (++f, a.doDasSearch(b.dasSource, e, t))
                })(this.tiers[q])
            }
        };
        n.prototype.doDasSearch = function(h, e, a) {
            h.features(null, {
                group: e,
                type: "transcript"
            }, function(b) {
                b || (b = []);
                for (var h = [], g = 0; g < b.length; ++g) {
                    var l = b[g];
                    l.label.toLowerCase() == e.toLowerCase() && h.push(l)
                }
                return a(h)
            }, !1)
        }
    }, {
        "./bin": 4,
        "./cbrowser": 6,
        "./trix": 46
    }],
    31: [function(e, u, s) {
        function p(h, g) {
            g || (g = a);
            for (var l = b.length; h * b[l % b.length] * Math.pow(10, l / b.length | 0) < g;) ++l;
            return b[l % b.length] * Math.pow(10, l / b.length | 0)
        }

        function n(a, b) {
            var h = a.viewport.getContext("2d"),
                f = a.browser.retina && (1 < window.devicePixelRatio || 1 < a.browser.devicePixelRatio),
                d = a.browser.featurePanelWidth + 2E3;
            f && (d *= 2);
            var e = a.viewport.width | 0;
            e < d - 50 && (a.viewport.width = e = d);
            d = 50;
            b && b.seq && (d += 25);
            var n = d;
            f && (n *= 2);
            a.viewport.height = n;
            a.viewport.style.height = "" + d + "px";
            a.viewport.style.width = f ? "" + e / 2 + "px" : "" + e + "px";
            a.layoutHeight = d;
            a.updateHeight();
            a.background && (h.fillStyle = a.background, h.fillRect(0, 0, e, a.viewport.height));
            f && h.scale(2, 2);
            h.translate(1E3, 0);
            m(a, b, h);
            a.norigin = a.browser.viewStart;
            a.viewportHolder.style.left = "-1000px"
        }

        function m(a, b, h) {
            var f = a.browser.scale,
                d = a.browser.viewStart - 1E3 / f | 0,
                e = a.browser.viewEnd + 2E3 / f,
                m = a.browser.currentSeqMax,
                n = e;
            0 < m && m < e && (n = m);
            for (var r = p(f), s = Math.max(0, (d / r | 0) * r), m = a.browser.viewStart; s <= n;) h.fillStyle = 0 == s / r % 2 ? "white" : "black", h.strokeStyle = "black", h.fillRect((s - m) * f, 8, r * f, 3), h.strokeRect((s - m) * f, 8, r * f, 3), h.fillStyle = "black", h.fillText(v(s), (s - m) * f, 22), s += r;
            if (b && b.seq)
                for (; d <= e; ++d) d >= b.start && d <= b.end && (n = b.seq.substr(d - b.start, 1).toUpperCase(), (r = a.browser.baseColors[n]) || (r = "gray"), h.fillStyle = r, 8 <= f ? (r = h.measureText(n).width, h.fillText(n, (d - m) * f + 0.5 * (f - r), 52)) : h.fillRect((d - m) * f, 42, f, 12))
        }

        function h(a, b) {
            var h = a.browser.scale,
                f = a.browser.viewStart - 1E3 / h | 0,
                d = a.browser.viewEnd + 2E3 / h,
                e = a.browser.currentSeqMax,
                m = d;
            0 < e && e < d && (m = e);
            for (var n = p(h), s = Math.max(0, (f / n | 0) * n), e = a.browser.viewStart, L = r(z, "g", [], {
                    fontSize: "8pt"
                }); s <= m;) L.appendChild(r(z, "rect", null, {
                x: (s - e) * h,
                y: 8,
                width: n * h,
                height: 3,
                fill: 0 == s / n % 2 ? "white" : "black",
                stroke: "black"
            })), L.appendChild(r(z,
                "text", v(s), {
                    x: (s - e) * h,
                    y: 28,
                    fill: "black",
                    stroke: "none"
                })), s += n;
            if (b && b.seq)
                for (; f <= d; ++f) f >= b.start && f <= b.end && (m = b.seq.substr(f - b.start, 1).toUpperCase(), (n = a.browser.baseColors[m]) || (n = "gray"), 8 <= h ? L.appendChild(r(z, "text", m, {
                    x: (0.5 + f - e) * h,
                    y: 52,
                    textAnchor: "middle",
                    fill: n
                })) : L.appendChild(r(z, "rect", null, {
                    x: (f - e) * h,
                    y: 42,
                    width: h,
                    height: 12,
                    fill: n
                })));
            return L
        }
        if ("undefined" !== typeof e) {
            s = e("./utils");
            var v = s.formatLongInt,
                r = s.makeElementNS,
                z = e("./svg-utils").NS_SVG,
                v = e("./numformats").formatLongInt
        }
        var a =
            100,
            b = [1, 2, 5],
            z = "http://www.w3.org/2000/svg";
        "undefined" !== typeof u && (u.exports = {
            drawSeqTier: n,
            drawSeqTierGC: m,
            svgSeqTier: h
        })
    }, {
        "./numformats": 26,
        "./svg-utils": 38,
        "./utils": 48
    }],
    32: [function(e, u, s) {
        if ("undefined" != typeof e) {
            u = e("./cbrowser");
            var p = u.Browser,
                n = u.sourceDataURI,
                m = u.sourcesAreEqual,
                h = e("./version"),
                v = e("./utils").miniJSONify,
                r = e("./sha1").hex_sha1
        }
        p.prototype.nukeStatus = function() {
            delete localStorage["dalliance." + this.cookieKey + ".view-chr"];
            delete localStorage["dalliance." + this.cookieKey +
                ".view-start"];
            delete localStorage["dalliance." + this.cookieKey + ".view-end"];
            delete localStorage["dalliance." + this.cookieKey + ".current-seq-length"];
            delete localStorage["dalliance." + this.cookieKey + ".showing-alt-zoom"];
            delete localStorage["dalliance." + this.cookieKey + ".saved-zoom"];
            delete localStorage["dalliance." + this.cookieKey + ".sources"];
            delete localStorage["dalliance." + this.cookieKey + ".hubs"];
            delete localStorage["dalliance." + this.cookieKey + ".version"];
            delete localStorage["dalliance." + this.cookieKey +
                ".reverse-scrolling"];
            delete localStorage["dalliance." + this.cookieKey + ".reverse-key-scrolling"];
            delete localStorage["dalliance." + this.cookieKey + ".ruler-location"]
        };
        p.prototype.storeStatus = function() {
            this.storeViewStatus();
            this.storeTierStatus()
        };
        p.prototype.storeViewStatus = function() {
            !this.cookieKey || this.noPersist || this.noPersistView || (localStorage["dalliance." + this.cookieKey + ".view-chr"] = this.chr, localStorage["dalliance." + this.cookieKey + ".view-start"] = this.viewStart | 0, localStorage["dalliance." +
                this.cookieKey + ".view-end"] = this.viewEnd | 0, localStorage["dalliance." + this.cookieKey + ".showing-alt-zoom"] = "" + this.isSnapZooming, localStorage["dalliance." + this.cookieKey + ".saved-zoom"] = this.savedZoom, this.currentSeqMax && (localStorage["dalliance." + this.cookieKey + ".current-seq-length"] = this.currentSeqMax))
        };
        p.prototype.storeTierStatus = function() {
            if (this.cookieKey && !this.noPersist) {
                for (var e = [], a = 0; a < this.tiers.length; ++a) {
                    var b = this.tiers[a];
                    b.dasSource.noPersist || e.push({
                        source: b.dasSource,
                        config: b.config || {}
                    })
                }
                localStorage["dalliance." + this.cookieKey + ".sources"] = JSON.stringify(e);
                e = {};
                a = [];
                for (b = 0; b < this.hubObjects.length; ++b) {
                    var m = this.hubObjects[b],
                        g = {
                            url: m.hub.url,
                            genome: m.genome
                        };
                    m.credentials && (g.credentials = m.credentials);
                    m.mapping && (g.mapping = m.mapping);
                    e[g.url] = !0;
                    a.push(g)
                }
                for (b = 0; b < this.hubs.length; ++b) g = this.hubs[b], "string" === typeof g && (g = {
                    url: g
                }), e[g.url] || a.push(g);
                localStorage["dalliance." + this.cookieKey + ".hubs"] = JSON.stringify(a);
                localStorage["dalliance." + this.cookieKey + ".reverse-scrolling"] =
                    this.reverseScrolling;
                localStorage["dalliance." + this.cookieKey + ".reverse-key-scrolling"] = this.reverseKeyScrolling;
                localStorage["dalliance." + this.cookieKey + ".single-base-highlight"] = this.singleBaseHighlight;
                localStorage["dalliance." + this.cookieKey + ".ruler-location"] = this.rulerLocation;
                localStorage["dalliance." + this.cookieKey + ".export-ruler"] = this.exportRuler;
                localStorage["dalliance." + this.cookieKey + ".export-highlights"] = this.exportHighlights;
                localStorage["dalliance." + this.cookieKey + ".version"] = h.CONFIG
            }
        };
        p.prototype.restoreStatus = function() {
            if (!this.noPersist) {
                var e = localStorage["dalliance." + this.cookieKey + ".version"];
                if (h.CONFIG == (e ? e | 0 : -100)) {
                    var e = localStorage["dalliance." + this.cookieKey + ".configHash"] || "",
                        a = r(v({
                            sources: this.sources,
                            hubs: this.hubs
                        }));
                    if (a != e) localStorage["dalliance." + this.cookieKey + ".configHash"] = a;
                    else {
                        e = {};
                        for (a = 0; a < this.sources.length; ++a) {
                            var b = this.sources[a];
                            if (b) {
                                var q = n(b),
                                    g = e[q];
                                g || (e[q] = g = []);
                                g.push(b)
                            }
                        }
                        if (!this.noPersistView && (a = localStorage["dalliance." + this.cookieKey +
                                ".view-chr"], b = localStorage["dalliance." + this.cookieKey + ".view-start"] | 0, g = localStorage["dalliance." + this.cookieKey + ".view-end"] | 0, a && b && g)) {
                            this.chr = a;
                            this.viewStart = b;
                            this.viewEnd = g;
                            if (a = localStorage["dalliance." + this.cookieKey + ".current-seq-length"]) this.currentSeqMax = a | 0;
                            this.isSnapZooming = "true" == localStorage["dalliance." + this.cookieKey + ".showing-alt-zoom"];
                            a = parseFloat(localStorage["dalliance." + this.cookieKey + ".saved-zoom"]);
                            "number" !== typeof a || isNaN(a) || (this.savedZoom = a)
                        }
                        this.reverseScrolling =
                            (a = localStorage["dalliance." + this.cookieKey + ".reverse-scrolling"]) && "true" == a;
                        this.reverseKeyScrolling = (a = localStorage["dalliance." + this.cookieKey + ".reverse-key-scrolling"]) && "true" == a;
                        this.singleBaseHighlight = (a = localStorage["dalliance." + this.cookieKey + ".single-base-highlight"]) && "true" == a;
                        if (a = localStorage["dalliance." + this.cookieKey + ".ruler-location"]) this.rulerLocation = a;
                        if (a = localStorage["dalliance." + this.cookieKey + ".export-ruler"]) this.exportRuler = "true" === a;
                        if (a = localStorage["dalliance." + this.cookieKey +
                                ".export-highlights"]) this.exportHighlights = "true" === a;
                        if (a = localStorage["dalliance." + this.cookieKey + ".sources"]) {
                            var l = JSON.parse(a);
                            this.sources = [];
                            this.restoredConfigs = [];
                            for (a = 0; a < l.length; ++a)
                                for (b = this.sources[a] = l[a].source, this.restoredConfigs[a] = l[a].config, q = n(b), g = e[q] || [], q = 0; q < g.length; ++q) {
                                    var f = g[q];
                                    if (m(b, f))
                                        for (var d in f) f.hasOwnProperty(d) && ("function" === typeof f[d] || f[d] instanceof Blob) && (b[d] = f[d])
                                }
                        }
                        if (d = localStorage["dalliance." + this.cookieKey + ".hubs"]) this.hubs = JSON.parse(d);
                        return !0
                    }
                }
            }
        };
        p.prototype.reset = function() {
            for (var h = this.tiers.length - 1; 0 <= h; --h) this.removeTier({
                index: h
            }, !0);
            for (h = 0; h < this.defaultSources.length; ++h) this.defaultSources[h].disabled || this.addTier(this.defaultSources[h]);
            this.highlights.splice(0, this.highlights.length);
            this.setLocation(this.defaultChr, this.defaultStart, this.defaultEnd)
        }
    }, {
        "./cbrowser": 6,
        "./sha1": 33,
        "./utils": 48,
        "./version": 50
    }],
    33: [function(e, u, s) {
        function p(a) {
            a = m(h(a));
            for (var b = r ? "0123456789ABCDEF" : "0123456789abcdef", e = "", g, l =
                    0; l < a.length; l++) g = a.charCodeAt(l), e += b.charAt(g >>> 4 & 15) + b.charAt(g & 15);
            return e
        }

        function n(a) {
            a = m(h(a));
            for (var b = "", e = a.length, g = 0; g < e; g += 3)
                for (var l = a.charCodeAt(g) << 16 | (g + 1 < e ? a.charCodeAt(g + 1) << 8 : 0) | (g + 2 < e ? a.charCodeAt(g + 2) : 0), f = 0; 4 > f; f++) b = 8 * g + 6 * f > 8 * a.length ? b + z : b + "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/".charAt(l >>> 6 * (3 - f) & 63);
            return b
        }

        function m(a) {
            for (var b = Array(a.length >> 2), h = 0; h < b.length; h++) b[h] = 0;
            for (h = 0; h < 8 * a.length; h += 8) b[h >> 5] |= (a.charCodeAt(h / 8) & 255) << 24 -
                h % 32;
            a = 8 * a.length;
            b[a >> 5] |= 128 << 24 - a % 32;
            b[(a + 64 >> 9 << 4) + 15] = a;
            a = Array(80);
            for (var h = 1732584193, g = -271733879, e = -1732584194, f = 271733878, d = -1009589776, m = 0; m < b.length; m += 16) {
                for (var n = h, p = g, r = e, s = f, u = d, z = 0; 80 > z; z++) {
                    if (16 > z) a[z] = b[m + z];
                    else {
                        var y = a[z - 3] ^ a[z - 8] ^ a[z - 14] ^ a[z - 16];
                        a[z] = y << 1 | y >>> 31
                    }
                    var y = h << 5 | h >>> 27,
                        H;
                    H = 20 > z ? g & e | ~g & f : 40 > z ? g ^ e ^ f : 60 > z ? g & e | g & f | e & f : g ^ e ^ f;
                    y = v(v(y, H), v(v(d, a[z]), 20 > z ? 1518500249 : 40 > z ? 1859775393 : 60 > z ? -1894007588 : -899497514));
                    d = f;
                    f = e;
                    e = g << 30 | g >>> 2;
                    g = h;
                    h = y
                }
                h = v(h, n);
                g = v(g, p);
                e = v(e, r);
                f = v(f, s);
                d = v(d, u)
            }
            b = [h, g, e, f, d];
            a = "";
            for (h = 0; h < 32 * b.length; h += 8) a += String.fromCharCode(b[h >> 5] >>> 24 - h % 32 & 255);
            return a
        }

        function h(a) {
            for (var b = "", h = -1, g, e; ++h < a.length;) g = a.charCodeAt(h), e = h + 1 < a.length ? a.charCodeAt(h + 1) : 0, 55296 <= g && 56319 >= g && 56320 <= e && 57343 >= e && (g = 65536 + ((g & 1023) << 10) + (e & 1023), h++), 127 >= g ? b += String.fromCharCode(g) : 2047 >= g ? b += String.fromCharCode(192 | g >>> 6 & 31, 128 | g & 63) : 65535 >= g ? b += String.fromCharCode(224 | g >>> 12 & 15, 128 | g >>> 6 & 63, 128 | g & 63) : 2097151 >= g && (b += String.fromCharCode(240 | g >>> 18 &
                7, 128 | g >>> 12 & 63, 128 | g >>> 6 & 63, 128 | g & 63));
            return b
        }

        function v(a, b) {
            var h = (a & 65535) + (b & 65535);
            return (a >> 16) + (b >> 16) + (h >> 16) << 16 | h & 65535
        }
        var r = 0,
            z = "";
        "undefined" !== typeof u && (u.exports = {
            b64_sha1: n,
            hex_sha1: p
        })
    }, {}],
    34: [function(e, u, s) {
        function p(a, c) {
            ka[a] = c
        }

        function n(a, c) {
            ga[a] = c
        }

        function m(a) {
            if (ga[a]) return ga[a](a)
        }

        function h(a) {
            var c = this;
            this.source = a;
            this.cfsid = "cfs" + ++la;
            a.name && (this.name = a.name);
            a.addChangeListener && a.addChangeListener(function() {
                c.cfsid = "cfs" + ++la
            })
        }

        function v() {
            this.busy =
                0;
            this.activityListeners = [];
            this.readinessListeners = [];
            this.readiness = null
        }

        function r(a) {
            this.dasSource = new w(a);
            this.busy = 0;
            this.activityListeners = []
        }

        function z(a) {
            this.dasSource = new w(a);
            this.awaitedEntryPoints = new N;
            var c = this;
            this.dasSource.entryPoints(function(a) {
                c.awaitedEntryPoints.provide(a)
            })
        }

        function a(a) {
            var c = this;
            this.source = a;
            this.twoBit = new N;
            if (a.twoBitURI) a = new A(a.twoBitURI);
            else if (a.twoBitBlob) a = new J(a.twoBitBlob);
            else throw Error("No twoBitURI or twoBitBlob parameter");
            k(a, function(a,
                b) {
                b ? console.log(b) : c.twoBit.provide(a)
            })
        }

        function b(a) {
            v.call(this);
            var c = this;
            this.readiness = "Connecting";
            this.bwgSource = this.opts = a;
            c.bwgHolder = new N;
            if (this.opts.preflight) {
                var b = ea[this.opts.preflight];
                if (!b) {
                    b = new N;
                    ea[this.opts.preflight] = b;
                    var d = new XMLHttpRequest;
                    d.onreadystatechange = function() {
                        4 == d.readyState && (200 == d.status ? b.provide("success") : b.provide("failure"))
                    };
                    d.open("get", this.opts.preflight + "?" + hex_sha1("salt" + Date.now()), !0);
                    this.opts.credentials && (d.withCredentials = !0);
                    d.send("")
                }
                b.await(function(a) {
                    "success" ===
                    a && c.init()
                })
            } else c.init()
        }

        function q(a, c) {
            v.call(this);
            this.worker = c;
            this.readiness = "Connecting";
            this.bwgSource = this.opts = a;
            this.keyHolder = new N;
            this.init()
        }

        function g(a, c) {
            if (!(a.flag & K.SEGMENT_UNMAPPED)) {
                var b;
                b = a.seq ? a.seq.length : a.seqLength;
                if (a.cigar) {
                    b = 0;
                    for (var d = Z(a.cigar), f = 0; f < d.length; ++f) {
                        var g = d[f];
                        if ("M" == g.op || "D" == g.op) b += g.cnt
                    }
                }
                d = new I;
                d.min = a.pos + 1;
                d.max = a.pos + b;
                d.segment = a.segment;
                d.type = "bam";
                d.id = a.readName;
                d.notes = ["MQ=" + a.mq];
                d.cigar = a.cigar;
                d.seq = a.seq;
                d.quals = a.quals;
                d.orientation =
                    a.flag & K.REVERSE_COMPLEMENT ? "-" : "+";
                d.bamRecord = a;
                c && a.flag & K.MULTIPLE_SEGMENTS && (d.groups = [{
                    id: a.readName,
                    type: "readpair"
                }]);
                return d
            }
        }

        function l(a) {
            v.call(this);
            var c = this;
            this.bamSource = a;
            this.opts = {
                credentials: a.credentials,
                preflight: a.preflight,
                bamGroup: a.bamGroup
            };
            this.bamHolder = new N;
            if (this.opts.preflight) {
                var b = ea[this.opts.preflight];
                if (!b) {
                    b = new N;
                    ea[this.opts.preflight] = b;
                    var d = new XMLHttpRequest;
                    d.onreadystatechange = function() {
                        4 == d.readyState && (200 == d.status ? b.provide("success") : b.provide("failure"))
                    };
                    d.open("get", this.opts.preflight + "?" + hex_sha1("salt" + Date.now()), !0);
                    this.opts.credentials && (d.withCredentials = "true");
                    d.send("")
                }
                b.await(function(a) {
                    "success" === a && c.init()
                })
            } else c.init()
        }

        function f(a, c) {
            v.call(this);
            this.bamSource = a;
            this.worker = c;
            this.opts = {
                credentials: a.credentials,
                preflight: a.preflight,
                bamGroup: a.bamGroup
            };
            this.keyHolder = new N;
            this.init()
        }

        function d(a, c) {
            this.source = a;
            this.mapping = c;
            this.activityListeners = [];
            this.busy = 0
        }

        function t() {}

        function D() {}

        function P(a) {
            this.store = new T(a.jbURI,
                a.jbQuery)
        }
        if ("undefined" !== typeof e) {
            var G = e("./cbrowser").Browser,
                L = e("./tier").DasTier;
            s = e("./utils");
            var N = s.Awaited,
                Q = s.arrayIndexOf,
                y = s.shallowCopy,
                H = s.resolveUrlToPage;
            s = e("./das");
            var F = s.DASStylesheet,
                c = s.DASStyle,
                w = s.DASSource,
                B = s.DASSegment,
                I = s.DASFeature,
                x = s.DASSequence,
                C = s.DASLink;
            s = e("./bin");
            var A = s.URLFetchable,
                J = s.BlobFetchable,
                k = e("./twoBit").makeTwoBit,
                E = e("./bigwig").makeBwg;
            s = e("./bam");
            var S = s.makeBam,
                K = s.BamFlags;
            s = e("./spans");
            var R = s.Range,
                V = s.union,
                Z = e("./cigar").parseCigar,
                ja = e("./overlay").OverlayFeatureSource,
                T = e("./jbjson").JBrowseStore,
                ca = e("./chainset").Chainset;
            e("./style");
            var X = e("./encode").EncodeFetchable
        }
        var ka = {},
            ga = {};
        L.prototype.initSources = function() {
            var a = this,
                c = this.browser.createSources(this.dasSource);
            this.featureSource = c.features || new t;
            this.sequenceSource = c.sequence;
            this.featureSource && this.featureSource.addChangeListener && this.featureSource.addChangeListener(function() {
                a.browser.refreshTier(a)
            })
        };
        G.prototype.createSources = function(c) {
            var g = this.sourceCache.get(c);
            if (g) return g;
            var e, k;
            if ("sequence" == c.tier_type || c.twoBitURI || c.twoBitBlob) k = c.twoBitURI || c.twoBitBlob ? new a(c) : new z(c);
            else if (c.tier_type && ka[c.tier_type]) k = (0, ka[c.tier_type])(c), e = k.features, k = k.sequence;
            else if (c.bwgURI || c.bwgBlob) e = (e = this.getWorker()) ? new q(c, e) : new b(c);
            else if (c.bamURI || c.bamBlob) e = (e = this.getWorker()) ? new f(c, e) : new l(c);
            else if (c.jbURI) e = new P(c);
            else if (c.uri || c.features_uri) e = new r(c);
            if (c.overlay) {
                g = [];
                e && g.push(new h(e));
                for (e = 0; e < c.overlay.length; ++e) {
                    var m = this.createSources(c.overlay[e]);
                    m && m.features && g.push(m.features)
                }
                e = new ja(g, c)
            }
            c.sequenceAliases && (e = new d(e, new ca({
                type: "alias",
                sequenceAliases: c.sequenceAliases
            })));
            c.mapping && (e = new d(e, this.chains[c.mapping]));
            c.name && e && !e.name && (e.name = c.name);
            null != e && (e = new h(e));
            if (null != e || null != k) g = {
                features: e,
                sequence: k
            }, this.sourceCache.put(c, g);
            return g
        };
        L.prototype.fetchStylesheet = function(a) {
            (this.dasSource.stylesheet_uri || !(this.dasSource.tier_type || this.dasSource.bwgURI || this.dasSource.bwgBlob || this.dasSource.bamURI || this.dasSource.bamBlob ||
                this.dasSource.twoBitURI || this.dasSource.twoBitBlob || this.dasSource.jbURI || this.dasSource.overlay) ? new r(this.dasSource) : this.getSource()).getStyleSheet(a)
        };
        var la = 0;
        h.prototype.addReadinessListener = function(a) {
            if (this.source.addReadinessListener) return this.source.addReadinessListener(a);
            a(null)
        };
        h.prototype.search = function(a, c) {
            if (this.source.search) return this.source.search(a, c)
        };
        h.prototype.getDefaultFIPs = function(a) {
            if (this.source.getDefaultFIPs) return this.source.getDefaultFIPs(a)
        };
        h.prototype.getStyleSheet =
            function(a) {
                this.source.getStyleSheet(a)
            };
        h.prototype.getScales = function() {
            return this.source.getScales()
        };
        h.prototype.addActivityListener = function(a) {
            this.source.addActivityListener && this.source.addActivityListener(a)
        };
        h.prototype.addChangeListener = function(a) {
            this.source.addChangeListener && this.source.addChangeListener(a)
        };
        h.prototype.findNextFeature = function(a, c, b, d) {
            this.source.findNextFeature(a, c, b, d)
        };
        h.prototype.quantFindNextFeature = function(a, c, b, d, f) {
            this.source.quantFindNextFeature(a, c,
                b, d, f)
        };
        h.prototype.capabilities = function() {
            return this.source.capabilities ? this.source.capabilities() : {}
        };
        h.prototype.fetch = function(a, c, b, d, f, g, h, e) {
            if (!g) throw Error("Fetch pool is null");
            var k = this,
                l = this.cfsid,
                m = g.awaitedFeatures[l];
            if (m && m.started) {
                if (m.styleFilters.doesNotContain(e)) {
                    k.source.fetch(a, c, b, d, f, g, h, e);
                    return
                }
            } else m ? m.styleFilters.addAll(e) : (m = new N, m.styleFilters = e, g.awaitedFeatures[l] = m, g.requestsIssued.then(function() {
                m.started = !0;
                k.source.fetch(a, c, b, d, m.styleFilters.typeList(),
                    g,
                    function(a, c, b, d) {
                        m.res || m.provide({
                            status: a,
                            features: c,
                            scale: b,
                            coverage: d
                        })
                    }, m.styleFilters)
            }).catch(function(a) {
                console.log(a)
            }));
            m.await(function(a) {
                h(a.status, a.features, a.scale, a.coverage)
            })
        };
        v.prototype.addReadinessListener = function(a) {
            this.readinessListeners.push(a);
            a(this.readiness)
        };
        v.prototype.notifyReadiness = function() {
            for (var a = 0; a < this.readinessListeners.length; ++a) try {
                this.readinessListeners[a](this.readiness)
            } catch (c) {
                console.log(c)
            }
        };
        v.prototype.addActivityListener = function(a) {
            this.activityListeners.push(a)
        };
        v.prototype.notifyActivity = function() {
            for (var a = 0; a < this.activityListeners.length; ++a) try {
                this.activityListeners[a](this.busy)
            } catch (c) {
                console.log(c)
            }
        };
        v.prototype.getScales = function() {
            return null
        };
        v.prototype.fetch = function(a, c, b, d, f, g, h) {
            return h(null, [], 1E9)
        };
        v.prototype.getStyleSheet = function(a) {
            var b = new F,
                d = new c;
            d.glyph = "BOX";
            d.BGCOLOR = "blue";
            d.FGCOLOR = "black";
            b.pushStyle({
                type: "default"
            }, null, d);
            return a(b)
        };
        r.prototype.addActivityListener = function(a) {
            this.activityListeners.push(a)
        };
        r.prototype.notifyActivity =
            function() {
                for (var a = 0; a < this.activityListeners.length; ++a) try {
                    this.activityListeners[a](this.busy)
                } catch (c) {
                    console.log(c)
                }
            };
        r.prototype.getStyleSheet = function(a) {
            this.dasSource.stylesheet(function(c) {
                a(c)
            }, function() {
                a(null, "Couldn't fetch DAS stylesheet")
            })
        };
        r.prototype.fetch = function(a, c, b, d, f, g, h) {
            if (f && 0 == f.length) h(null, [], d);
            else if (this.dasSource.uri || this.dasSource.features_uri) {
                if (this.dasSource.dasStaticFeatures && this.cachedStaticFeatures) return h(null, this.cachedStaticFeatures, this.cachedStaticScale);
                var e = !1 !== this.dasSource.maxbins;
                f = {
                    type: f
                };
                e && (f.maxbins = 1 + ((b - c) / d | 0));
                var k = this;
                k.busy++;
                k.notifyActivity();
                this.dasSource.features(new B(a, c, b), f, function(a, c) {
                    k.busy--;
                    k.notifyActivity();
                    var b = d;
                    e || (b = 0.1);
                    !c && k.dasSource.dasStaticFeatures && (k.cachedStaticFeatures = a, k.cachedStaticScale = b);
                    h(c, a, b)
                })
            }
        };
        r.prototype.findNextFeature = this.sourceFindNextFeature = function(a, c, b, d) {
            if (this.dasSource.capabilities && 0 <= Q(this.dasSource.capabilities, "das1:adjacent-feature")) {
                var f = this;
                if (this.dasAdjLock) return console.log("Already looking for a next feature, be patient!");
                this.dasAdjLock = !0;
                a = {
                    adjacent: a + ":" + (c | 0) + ":" + (0 < b ? "F" : "B")
                };
                if (c = thisTier.getDesiredTypes(thisTier.browser.scale)) a.types = c;
                thisTier.dasSource.features(null, a, function(a) {
                    f.dasAdjLock = !1;
                    0 < a.length && null != a[0] && d(a[0])
                })
            }
        };
        z.prototype.fetch = function(a, c, b, d, f) {
            this.dasSource.sequence(new B(a, c, b), function(a) {
                return 1 == a.length ? f(null, a[0]) : f("Didn't get sequence")
            })
        };
        z.prototype.getSeqInfo = function(a, c) {
            this.awaitedEntryPoints.await(function(b) {
                for (var d = 0; d < b.length; ++d)
                    if (b[d].name == a) return c({
                        length: b[d].end
                    });
                return c()
            })
        };
        a.prototype.fetch = function(a, c, b, d, f) {
            this.twoBit.await(function(d) {
                d.fetch(a, c, b, function(d, g) {
                    if (g) return f(g, null);
                    var h = new x(a, c, b, "DNA", d);
                    return f(null, h)
                })
            })
        };
        a.prototype.getSeqInfo = function(a, c) {
            this.twoBit.await(function(b) {
                b.getSeq(a) ? b.getSeq(a).length(function(a) {
                    c({
                        length: a
                    })
                }) : c()
            })
        };
        r.prototype.getScales = function() {
            return []
        };
        var ea = {};
        b.prototype = Object.create(v.prototype);
        b.prototype.init = function() {
            var a = this,
                c;
            c = (c = this.bwgSource.uri || this.bwgSource.bwgURI) ? "encode" ===
                this.bwgSource.transport ? new X(c, {
                    credentials: this.opts.credentials
                }) : new A(c, {
                    credentials: this.opts.credentials
                }) : new J(this.bwgSource.bwgBlob);
            E(c, function(c, b) {
                b ? (a.error = b, a.readiness = null, a.notifyReadiness(), a.bwgHolder.provide(null)) : (a.bwgHolder.provide(c), a.readiness = null, a.notifyReadiness(), "bigbed" == c.type && c.getExtraIndices(function(c) {
                    a.extraIndices = c
                }))
            })
        };
        b.prototype.capabilities = function() {
            var a = {
                leap: !0
            };
            this.bwgHolder.res && "bigwig" == this.bwgHolder.res.type && (a.quantLeap = !0);
            if (this.extraIndices &&
                0 < this.extraIndices.length) {
                a.search = [];
                for (var c = 0; c < this.extraIndices.length; ++c) a.search.push(this.extraIndices[c].field)
            }
            return a
        };
        b.prototype.fetch = function(a, c, b, d, f, g, h) {
            var e = this;
            this.bwgHolder.await(function(g) {
                if (null == g) return h(e.error || "Can't access binary file", null, null);
                var k;
                k = !f || 0 == f.length || 0 <= Q(f, "density");
                e.opts.clientBin && (k = !1);
                if ("bigwig" == g.type || k || "undefined" !== typeof e.opts.forceReduction) {
                    k = -1;
                    for (var l = 0; l < g.zoomLevels.length; ++l)
                        if (g.zoomLevels[l].reduction <= d) k =
                            l;
                        else break;
                        "undefined" !== typeof e.opts.forceReduction && (k = e.opts.forceReduction);
                    k = 0 > k ? g.getUnzoomedView() : g.getZoomedView(k)
                } else k = g.getUnzoomedView();
                e.busy++;
                e.notifyActivity();
                k.readWigData(a, c, b, function(a) {
                    e.busy--;
                    e.notifyActivity();
                    var d = 1E9;
                    if ("bigwig" === g.type) {
                        var f = (b - c) / a.length / 2;
                        f < d && (d = f)
                    }
                    if (e.opts.link)
                        for (f = 0; f < a.length; ++f) {
                            var k = a[f];
                            k.label && (k.links = [new C("Link", e.opts.link.replace(/\$\$/, k.label))])
                        }
                    h(null, a, d)
                })
            })
        };
        b.prototype.quantFindNextFeature = function(a, c, b, d, f) {
            var g =
                this;
            g.busy++;
            g.notifyActivity();
            this.bwgHolder.res.thresholdSearch(a, c, b, d, function(a, c) {
                g.busy--;
                g.notifyActivity();
                return f(a, c)
            })
        };
        b.prototype.findNextFeature = function(a, c, b, d) {
            var f = this;
            f.busy++;
            f.notifyActivity();
            this.bwgHolder.res.getUnzoomedView().getFirstAdjacent(a, c, b, function(a) {
                f.busy--;
                f.notifyActivity();
                0 < a.length && null != a[0] && d(a[0])
            })
        };
        b.prototype.getScales = function() {
            var a = this.bwgHolder.res;
            if (a) {
                for (var c = [1], b = 0; b < a.zoomLevels.length; ++b) c.push(a.zoomLevels[b].reduction);
                return c
            }
            return null
        };
        b.prototype.search = function(a, c) {
            return this.extraIndices && 0 != this.extraIndices.length ? this.extraIndices[0].lookup(a, c) : c(null, "No indices available")
        };
        b.prototype.getDefaultFIPs = function(a) {
            if (this.opts.noExtraFeatureInfo) return !0;
            this.bwgHolder.await(function(c) {
                c && c.schema && c.definedFieldCount < c.schema.fields.length && a(function(a, b) {
                    for (var d = 0; d < b.hit.length; ++d)
                        if (b.hit[d].isSuperGroup) return;
                    for (d = c.definedFieldCount; d < c.schema.fields.length; ++d) {
                        var f = c.schema.fields[d];
                        b.add(f.comment, a[f.name])
                    }
                })
            })
        };
        b.prototype.getStyleSheet = function(a) {
            this.bwgHolder.await(function(b) {
                if (!b) return a(null, "bbi error");
                var d = new F;
                if ("bigbed" == b.type) {
                    var f = new c;
                    f.glyph = "BOX";
                    f.FGCOLOR = "black";
                    f.BGCOLOR = "blue";
                    f.HEIGHT = 8;
                    f.BUMP = !0;
                    f.LABEL = !0;
                    f.ZINDEX = 20;
                    d.pushStyle({
                        type: "bigwig"
                    }, null, f);
                    f.glyph = "BOX";
                    f.FGCOLOR = "black";
                    f.BGCOLOR = "red";
                    f.HEIGHT = 10;
                    f.BUMP = !0;
                    f.ZINDEX = 20;
                    d.pushStyle({
                        type: "translation"
                    }, null, f);
                    f = new c;
                    f.glyph = "BOX";
                    f.FGCOLOR = "black";
                    f.BGCOLOR = "white";
                    f.HEIGHT = 10;
                    f.ZINDEX = 10;
                    f.BUMP = !0;
                    f.LABEL = !0;
                    d.pushStyle({
                        type: "transcript"
                    }, null, f);
                    f = new c;
                    f.glyph = "HISTOGRAM";
                    f.COLOR1 = "white";
                    f.COLOR2 = "black";
                    f.HEIGHT = 30;
                    d.pushStyle({
                        type: "density"
                    }, null, f)
                } else f = new c, f.glyph = "HISTOGRAM", f.COLOR1 = "white", f.COLOR2 = "black", f.HEIGHT = 30, d.pushStyle({
                    type: "default"
                }, null, f);
                12 == b.definedFieldCount && 14 <= b.fieldCount && (d.geneHint = !0);
                return a(d)
            })
        };
        q.prototype = Object.create(v.prototype);
        q.prototype.init = function() {
            var a = this,
                c = this.uri || this.bwgSource.uri || this.bwgSource.bwgURI,
                b = this.bwgSource.blob || this.bwgSource.bwgBlob,
                d = function(c, b) {
                    a.readiness = null;
                    a.notifyReadiness();
                    c ? a.worker.postCommand({
                        command: "meta",
                        connection: c
                    }, function(b, d) {
                        d ? (a.error = d, a.keyHolder.provide(null)) : (a.meta = b, a.keyHolder.provide(c))
                    }) : (a.error = b, a.keyHolder.provide(null))
                };
            b ? this.worker.postCommand({
                command: "connectBBI",
                blob: b
            }, d) : this.worker.postCommand({
                command: "connectBBI",
                uri: H(c),
                transport: this.bwgSource.transport,
                credentials: this.bwgSource.credentials
            }, d)
        };
        q.prototype.capabilities = function() {
            var a = {
                leap: !0
            };
            this.meta && "bigwig" ==
                this.meta.type && (a.quantLeap = !0);
            if (this.meta && this.meta.extraIndices && 0 < this.meta.extraIndices.length) {
                a.search = [];
                for (var c = 0; c < this.meta.extraIndices.length; ++c) a.search.push(this.meta.extraIndices[c].field)
            }
            return a
        };
        q.prototype.fetch = function(a, c, b, d, f, g, h) {
            var e = this;
            e.busy++;
            e.notifyActivity();
            this.keyHolder.await(function(g) {
                if (!g) return e.busy--, e.notifyActivity(), h(e.error || "Can't access binary file", null, null);
                var k = -1,
                    l = !f || 0 == f.length || 0 <= Q(f, "density");
                e.opts.clientBin && (l = !1);
                if ("bigwig" ==
                    e.meta.type || l || "undefined" !== typeof e.opts.forceReduction) {
                    for (l = 1; l < e.meta.zoomLevels.length; ++l)
                        if (e.meta.zoomLevels[l] <= d) k = l - 1;
                        else break;
                        "undefined" !== typeof e.opts.forceReduction && (k = e.opts.forceReduction)
                }
                e.worker.postCommand({
                    command: "fetch",
                    connection: g,
                    chr: a,
                    min: c,
                    max: b,
                    zoom: k
                }, function(a, d) {
                    e.busy--;
                    e.notifyActivity();
                    var f = 1E9;
                    if ("bigwig" === e.meta.type) {
                        var g = (b - c) / a.length / 2;
                        g < f && (f = g)
                    }
                    if (e.opts.link)
                        for (g = 0; g < a.length; ++g) {
                            var k = a[g];
                            k.label && (k.links = [new C("Link", e.opts.link.replace(/\$\$/,
                                k.label))])
                        }
                    h(d, a, f)
                })
            })
        };
        q.prototype.quantFindNextFeature = function(a, c, b, d, f) {
            var g = this;
            this.busy++;
            this.notifyActivity();
            this.worker.postCommand({
                command: "quantLeap",
                connection: this.keyHolder.res,
                chr: a,
                pos: c,
                dir: b,
                threshold: d,
                under: !1
            }, function(a, c) {
                console.log(a, c);
                g.busy--;
                g.notifyActivity();
                return f(a, c)
            })
        };
        q.prototype.findNextFeature = function(a, c, b, d) {
            var f = this;
            this.busy++;
            this.notifyActivity();
            this.worker.postCommand({
                command: "leap",
                connection: this.keyHolder.res,
                chr: a,
                pos: c,
                dir: b
            }, function(a,
                c) {
                f.busy--;
                f.notifyActivity();
                0 < a.length && null != a[0] && d(a[0])
            })
        };
        q.prototype.getScales = function() {
            var a = this.meta;
            return a ? a.zoomLevels : null
        };
        q.prototype.search = function(a, c) {
            if (!this.meta.extraIndices || 0 == this.meta.extraIndices.length) return c(null, "No indices available");
            var b = this;
            this.busy++;
            this.notifyActivity();
            this.worker.postCommand({
                command: "search",
                connection: this.keyHolder.res,
                query: a,
                index: this.meta.extraIndices[0]
            }, function(a, d) {
                b.busy--;
                b.notifyActivity();
                c(a, d)
            })
        };
        q.prototype.getDefaultFIPs =
            function(a) {
                if (this.opts.noExtraFeatureInfo) return !0;
                var c = this;
                this.keyHolder.await(function(b) {
                    var d = c.meta;
                    d && d.schema && d.definedFieldCount < d.schema.fields.length && a(function(a, c) {
                        for (var b = 0; b < c.hit.length; ++b)
                            if (c.hit[b].isSuperGroup) return;
                        for (b = d.definedFieldCount; b < d.schema.fields.length; ++b) {
                            var f = d.schema.fields[b];
                            c.add(f.comment, a[f.name])
                        }
                    })
                })
            };
        q.prototype.getStyleSheet = function(a) {
            var b = this;
            this.keyHolder.await(function(d) {
                d = b.meta;
                if (!d) return a(null, "bbi error");
                var f = new F;
                if ("bigbed" ==
                    d.type) {
                    var g = new c;
                    g.glyph = "BOX";
                    g.FGCOLOR = "black";
                    g.BGCOLOR = "blue";
                    g.HEIGHT = 8;
                    g.BUMP = !0;
                    g.LABEL = !0;
                    g.ZINDEX = 20;
                    f.pushStyle({
                        type: "bigwig"
                    }, null, g);
                    g.glyph = "BOX";
                    g.FGCOLOR = "black";
                    g.BGCOLOR = "red";
                    g.HEIGHT = 10;
                    g.BUMP = !0;
                    g.ZINDEX = 20;
                    f.pushStyle({
                        type: "translation"
                    }, null, g);
                    g = new c;
                    g.glyph = "BOX";
                    g.FGCOLOR = "black";
                    g.BGCOLOR = "white";
                    g.HEIGHT = 10;
                    g.ZINDEX = 10;
                    g.BUMP = !0;
                    g.LABEL = !0;
                    f.pushStyle({
                        type: "transcript"
                    }, null, g);
                    g = new c;
                    g.glyph = "HISTOGRAM";
                    g.COLOR1 = "white";
                    g.COLOR2 = "black";
                    g.HEIGHT = 30;
                    f.pushStyle({
                            type: "density"
                        },
                        null, g)
                } else g = new c, g.glyph = "HISTOGRAM", g.COLOR1 = "white", g.COLOR2 = "black", g.HEIGHT = 30, f.pushStyle({
                    type: "default"
                }, null, g);
                12 == d.definedFieldCount && 14 <= d.fieldCount && (f.geneHint = !0);
                return a(f)
            })
        };
        l.prototype = Object.create(v.prototype);
        l.prototype.init = function() {
            var a = this,
                c, b;
            this.bamSource.bamBlob ? (c = new J(this.bamSource.bamBlob), b = new J(this.bamSource.baiBlob)) : (c = new A(this.bamSource.bamURI, {
                credentials: this.opts.credentials
            }), b = new A(this.bamSource.baiURI || this.bamSource.bamURI + ".bai", {
                credentials: this.opts.credentials
            }));
            S(c, b, null, function(c, b) {
                a.readiness = null;
                a.notifyReadiness();
                c ? a.bamHolder.provide(c) : (a.error = b, a.bamHolder.provide(null))
            })
        };
        l.prototype.fetch = function(a, c, b, d, f, h, e) {
            var k = f && 1 == f.length && "density" == f[0],
                l = this;
            l.busy++;
            l.notifyActivity();
            this.bamHolder.await(function(d) {
                if (!d) return l.busy--, l.notifyActivity(), e(l.error || "Couldn't fetch BAM");
                d.fetch(a, c, b, function(a, c) {
                    l.busy--;
                    l.notifyActivity();
                    if (c) e(c, null, null);
                    else {
                        for (var b = [], d = 0; d < a.length; ++d) {
                            var f = g(a[d], l.opts.bamGroup);
                            f && b.push(f)
                        }
                        e(null,
                            b, 1E9)
                    }
                }, {
                    light: k
                })
            })
        };
        l.prototype.getScales = function() {
            return 1E9
        };
        l.prototype.getStyleSheet = function(a) {
            this.bamHolder.await(function(b) {
                b = new F;
                var d = new c;
                d.glyph = "HISTOGRAM";
                d.COLOR1 = "black";
                d.COLOR2 = "red";
                d.HEIGHT = 30;
                b.pushStyle({
                    type: "density"
                }, "low", d);
                b.pushStyle({
                    type: "density"
                }, "medium", d);
                d = new c;
                d.glyph = "__SEQUENCE";
                d.FGCOLOR = "black";
                d.BGCOLOR = "blue";
                d.HEIGHT = 8;
                d.BUMP = !0;
                d.LABEL = !1;
                d.ZINDEX = 20;
                b.pushStyle({
                    type: "bam"
                }, "high", d);
                return a(b)
            })
        };
        f.prototype = Object.create(v.prototype);
        f.prototype.init =
            function() {
                var a = this,
                    c = this.bamSource.uri || this.bamSource.bamURI,
                    b = this.bamSource.indexUri || this.bamSource.baiURI || c + ".bai",
                    d = this.bamSource.bamBlob || this.bamSource.blob,
                    f = this.bamSource.baiBlob || this.bamSource.indexBlob,
                    g = function(c, b) {
                        a.readiness = null;
                        a.notifyReadiness();
                        c ? a.keyHolder.provide(c) : (a.error = b, a.keyHolder.provide(null))
                    };
                d ? this.worker.postCommand({
                    command: "connectBAM",
                    blob: d,
                    indexBlob: f
                }, g) : this.worker.postCommand({
                    command: "connectBAM",
                    uri: H(c),
                    indexUri: H(b),
                    credentials: this.bamSource.credentials,
                    indexChunks: this.bamSource.indexChunks
                }, g)
            };
        f.prototype.fetch = function(a, c, b, d, f, h, e) {
            var k = f && 1 == f.length && "density" == f[0],
                l = this;
            l.busy++;
            l.notifyActivity();
            this.keyHolder.await(function(d) {
                if (!d) return l.busy--, l.notifyActivity(), e(l.error || "Couldn't fetch BAM");
                l.worker.postCommand({
                    command: "fetch",
                    connection: d,
                    chr: a,
                    min: c,
                    max: b,
                    opts: {
                        light: k
                    }
                }, function(a, c) {
                    l.busy--;
                    l.notifyActivity();
                    if (c) e(c, null, null);
                    else {
                        for (var b = [], d = 0; d < a.length; ++d) {
                            var f = g(a[d], l.opts.bamGroup);
                            f && b.push(f)
                        }
                        e(null,
                            b, 1E9)
                    }
                })
            })
        };
        f.prototype.getScales = function() {
            return 1E9
        };
        f.prototype.getStyleSheet = function(a) {
            this.keyHolder.await(function(b) {
                b = new F;
                var d = new c;
                d.glyph = "HISTOGRAM";
                d.COLOR1 = "black";
                d.COLOR2 = "red";
                d.HEIGHT = 30;
                b.pushStyle({
                    type: "density"
                }, "low", d);
                b.pushStyle({
                    type: "density"
                }, "medium", d);
                d = new c;
                d.glyph = "__SEQUENCE";
                d.FGCOLOR = "black";
                d.BGCOLOR = "blue";
                d.HEIGHT = 8;
                d.BUMP = !0;
                d.LABEL = !1;
                d.ZINDEX = 20;
                b.pushStyle({
                    type: "bam"
                }, "high", d);
                return a(b)
            })
        };
        d.prototype.addActivityListener = function(a) {
            this.activityListeners.push(a)
        };
        d.prototype.notifyActivity = function() {
            for (var a = 0; a < this.activityListeners.length; ++a) try {
                this.activityListeners[a](this.busy)
            } catch (c) {
                console.log(c)
            }
        };
        d.prototype.getStyleSheet = function(a) {
            return this.source.getStyleSheet(a)
        };
        d.prototype.getScales = function() {
            return this.source.getScales()
        };
        d.prototype.getDefaultFIPs = function(a) {
            if (this.source.getDefaultFIPs) return this.source.getDefaultFIPs(a)
        };
        d.prototype.simplifySegments = function(a, c) {
            if (0 == a.length) return a;
            a.sort(function(a, c) {
                var b = a.name -
                    c.name;
                return b ? b : (b = a.start - c.start) ? b : a.end - c.end
            });
            for (var b = [], d = a[0], f = 0; f < a.length; ++f) {
                var g = a[f];
                g.name != d.name || g.start > d.end + c ? (b.push(d), d = g) : d = new B(d.name, Math.min(d.start, g.start), Math.max(d.end, g.end))
            }
            b.push(d);
            return b
        };
        d.prototype.fetch = function(a, c, b, d, f, g, h, e) {
            var k = this,
                l = b - c + 1;
            k.busy++;
            k.notifyActivity();
            this.mapping.sourceBlocksForRange(a, c, b, function(a) {
                if (0 == a.length) k.busy--, k.notifyActivity(), h("No mapping available for this regions", [], d);
                else {
                    a = k.simplifySegments(a, Math.max(100,
                        0.05 * l));
                    var c = [],
                        b = null,
                        m = a.length,
                        n;
                    a.map(function(a) {
                        k.source.fetch(a.name, a.start, a.end, d, f, g, function(d, f, g) {
                            d && !n && (n = d);
                            if (f)
                                for (d = 0; d < f.length; ++d) {
                                    var e = f[d],
                                        l = e.segment;
                                    0 == l.indexOf("chr") && (l = l.substr(3));
                                    l = k.mapping.mapSegment(l, e.min, e.max);
                                    if (0 == l.length) e.parts && 0 < e.parts.length && c.push(e);
                                    else
                                        for (var q = 0; q < l.length; ++q) {
                                            var w = l[q],
                                                p = y(e);
                                            p.segment = w.segment;
                                            p.min = w.min;
                                            p.max = w.max;
                                            w.partialMin && (p.partialMin = w.partialMin);
                                            w.partialMax && (p.partialMax = w.partialMax);
                                            w.flipped && ("-" ==
                                                e.orientation ? p.orientation = "+" : "+" == e.orientation && (p.orientation = "-"));
                                            c.push(p)
                                        }
                                }
                            f = k.mapping.mapPoint(a.name, a.start);
                            d = k.mapping.mapPoint(a.name, a.end);
                            f && d && (f = new R(f.pos, d.pos), b = b ? V(b, f) : f);
                            --m;
                            0 == m && (k.busy--, k.notifyActivity(), h(n, c, g, b))
                        }, e)
                    })
                }
            })
        };
        t.prototype.getScales = function() {
            return null
        };
        t.prototype.fetch = function(a, c, b, d, f, g, h) {
            return h(null, [], 1E9)
        };
        t.prototype.getStyleSheet = function(a) {
            var b = new F,
                d = new c;
            d.glyph = "BOX";
            d.BGCOLOR = "blue";
            d.FGCOLOR = "black";
            b.pushStyle({
                    type: "default"
                },
                null, d);
            return a(b)
        };
        D.prototype.fetch = function(a, c, b, d, f) {
            return f(null, null)
        };
        P.prototype.getScales = function() {
            return null
        };
        P.prototype.getStyleSheet = function(a) {
            var b = new F,
                d = new c;
            d.glyph = "BOX";
            d.FGCOLOR = "black";
            d.BGCOLOR = "green";
            d.HEIGHT = 8;
            d.BUMP = !0;
            d.LABEL = !0;
            d.ZINDEX = 20;
            b.pushStyle({
                type: "default"
            }, null, d);
            return a(b)
        };
        P.prototype.fetch = function(a, c, b, d, f, g, h) {
            f && 0 == f.length ? h(null, [], d) : this.store.features(new B(a, c, b), {}, function(a, c) {
                h(c, a, 1E5)
            })
        };
        G.prototype.sourceAdapterIsCapable = function(a,
            c) {
            return a.capabilities ? a.capabilities()[c] : !1
        };
        "undefined" !== typeof u && (u.exports = {
            FeatureSourceBase: v,
            TwoBitSequenceSource: a,
            DASSequenceSource: z,
            MappedFeatureSource: d,
            CachingFeatureSource: h,
            BWGFeatureSource: b,
            RemoteBWGFeatureSource: q,
            BAMFeatureSource: l,
            RemoteBAMFeatureSource: f,
            DummyFeatureSource: t,
            DummySequenceSource: D,
            registerSourceAdapterFactory: p,
            registerParserFactory: n,
            makeParser: m
        }, e("./ensembljson"), e("./tabix-source"), e("./memstore"), e("./bedwig"), e("./vcf"))
    }, {
        "./bam": 1,
        "./bedwig": 2,
        "./bigwig": 3,
        "./bin": 4,
        "./cbrowser": 6,
        "./chainset": 7,
        "./cigar": 8,
        "./das": 10,
        "./encode": 12,
        "./ensembljson": 13,
        "./jbjson": 22,
        "./memstore": 25,
        "./overlay": 27,
        "./spans": 35,
        "./style": 36,
        "./tabix-source": 39,
        "./tier": 44,
        "./twoBit": 47,
        "./utils": 48,
        "./vcf": 49
    }],
    35: [function(e, u, s) {
        function p(a, b) {
            if ("number" != typeof a || "number" != typeof b) throw "Bad range " + a + "," + b;
            this._min = a;
            this._max = b
        }

        function n(a) {
            this._ranges = a
        }

        function m(a, b) {
            a instanceof Array || (a = [a], b && a.push(b));
            if (0 == a.length) return null;
            if (1 == a.length) return a[0];
            for (var h = [], g = 0; g < a.length; ++g) a[g]._pushRanges(h);
            for (var h = h.sort(z), g = [], e = h[0], e = new p(e._min, e._max), f = 1; f < h.length; ++f) {
                var d = h[f];
                d._min > e._max + 1 ? (g.push(e), e = new p(d._min, d._max)) : d._max > e._max && (e._max = d._max)
            }
            g.push(e);
            return 1 == g.length ? g[0] : new n(g)
        }

        function h(a, b) {
            for (var h = a.ranges(), g = b.ranges(), e = h.length, f = g.length, d = 0, m = 0, r = []; d < e && m < f;) {
                a = h[d];
                b = g[m];
                var s = Math.max(a.min(), b.min()),
                    v = Math.min(a.max(), b.max());
                v >= s && r.push(new p(s, v));
                a.max() > b.max() ? ++m : ++d
            }
            return 0 == r.length ?
                null : 1 == r.length ? r[0] : new n(r)
        }

        function v(a) {
            var b = 0;
            a = a.ranges();
            for (var h = 0; h < a.length; ++h) var g = a[h],
                b = b + (g.max() - g.min() + 1);
            return b
        }

        function r(a, b) {
            return a.min() < b.min() ? -1 : a.min() > b.min() ? 1 : a.max() < b.max() ? -1 : b.max() > a.max() ? 1 : 0
        }

        function z(a, b) {
            return a._min < b._min ? -1 : a._min > b._min ? 1 : a._max < b._max ? -1 : b._max > a._max ? 1 : 0
        }
        p.prototype.min = function() {
            return this._min
        };
        p.prototype.max = function() {
            return this._max
        };
        p.prototype.contains = function(a) {
            return a >= this._min && a <= this._max
        };
        p.prototype.isContiguous =
            function() {
                return !0
            };
        p.prototype.ranges = function() {
            return [this]
        };
        p.prototype._pushRanges = function(a) {
            a.push(this)
        };
        p.prototype.toString = function() {
            return "[" + this._min + "-" + this._max + "]"
        };
        n.prototype.min = function() {
            return this._ranges[0].min()
        };
        n.prototype.max = function() {
            return this._ranges[this._ranges.length - 1].max()
        };
        n.prototype.contains = function(a) {
            for (var b = 0; b < this._ranges.length; ++b)
                if (this._ranges[b].contains(a)) return !0;
            return !1
        };
        n.prototype.isContiguous = function() {
            return 1 < this._ranges.length
        };
        n.prototype.ranges = function() {
            return this._ranges
        };
        n.prototype._pushRanges = function(a) {
            for (var b = 0; b < this._ranges.length; ++b) a.push(this._ranges[b])
        };
        n.prototype.toString = function() {
            for (var a = "", b = 0; b < this._ranges.length; ++b) 0 < b && (a += ","), a += this._ranges[b].toString();
            return a
        };
        "undefined" !== typeof u && (u.exports = {
            Range: p,
            union: m,
            intersection: h,
            coverage: v,
            rangeOver: r,
            _rangeOrder: z
        })
    }, {}],
    36: [function(e, u, s) {
        function p(e, h, n) {
            this.type = e;
            this.method = h;
            this.label = n
        }

        function n(e) {
            this._filters = {};
            if (e)
                for (var h =
                        0; h < e.length; ++h) this.add(e[h])
        }
        p.prototype.equals = function(e) {
            return this.type == e.type && this.method == e.method && this.label == e.label
        };
        p.prototype.toString = function() {
            var e = [];
            this.type && e.push("type=" + this.type);
            this.method && e.push("method=" + this.method);
            this.label && e.push("label=" + this.label);
            return "StyleFilter<" + e.join(";") + ">"
        };
        n.prototype.add = function(e) {
            var h = e.toString();
            this._filters[h] || (this._filters[h] = e, this._list = null)
        };
        n.prototype.addAll = function(e) {
            e = e.list();
            for (var h = 0; h < e.length; ++h) this.add(e[h])
        };
        n.prototype.doesNotContain = function(e) {
            e = e.list();
            for (var h = 0; h < e.length; ++h)
                if (!this._filters[h.toString()]) return !0;
            return !1
        };
        n.prototype.list = function() {
            if (!this._list) {
                this._list = [];
                for (var e in this._filters) this._filters.hasOwnProperty(e) && this._list.push(this._filters[e])
            }
            return this._list
        };
        n.prototype.typeList = function() {
            for (var e = [], h = this.list(), n = 0; n < h.length; ++n) {
                var p = h[n].type;
                if (!p || "default" == p) return null;
                0 > e.indexOf(p) && e.push(p)
            }
            return e
        };
        "undefined" !== typeof u && (u.exports = {
            StyleFilter: p,
            StyleFilterSet: n
        })
    }, {}],
    37: [function(e, u, s) {
        if ("undefined" !== typeof e) {
            var p = e("./cbrowser").Browser,
                n = e("./utils").makeElementNS,
                m = e("./version"),
                h = e("./sequence-draw").svgSeqTier;
            u = e("./svg-utils");
            var v = u.NS_SVG,
                r = u.NS_XLINK,
                z = u.SVGPath,
                a = e("./numformats").formatQuantLabel
        }
        p.prototype.makeSVG = function(b) {
            b = b || {};
            var e = b.minTierHeight || 20,
                g = document.implementation.createDocument(v, "svg", null),
                l = n(v, "g", null, {
                    fontFamily: "helvetica",
                    fontSize: "8pt"
                });
            g.documentElement.appendChild(l);
            var f = n(v, "a",
                //n(v, "text", "Graphics from Dalliance " + m, {
                n(v, "text", "", {
                    x: (this.featurePanelWidth + 200 + 20) / 2,
                    y: 30,
                    strokeWidth: 0,
                    fontSize: "12pt",
                    textAnchor: "middle",
                    fill: "blue"
                }));
            f.setAttribute("xmlns:xlink", r);
            f.setAttribute("xlink:href", "http://www.biodalliance.org/");
            l.appendChild(f);
            f = n(v, "rect", null, {
                x: 200,
                y: 50,
                width: this.featurePanelWidth,
                height: 1E5
            });
            f = n(v, "clipPath", f, {
                id: "featureClip"
            });
            l.appendChild(f);
            for (var f = 70, d = n(v, "g", null, {}), p = 0; p < this.tiers.length; ++p) {
                var s = this.tiers[p],
                    u = n(v, "g", null, {
                        clipPath: "url(#featureClip)",
                        clipRule: "nonzero"
                    }),
                    G = n(v, "g"),
                    L = f,
                    N = n(v, "rect", null, {
                        x: 0,
                        y: L,
                        width: "10000",
                        height: 50,
                        fill: s.background
                    });
                u.appendChild(N);
                if (s.sequenceSource) {
                    var Q = h(s, s.currentSequence);
                    u.appendChild(n(v, "g", Q, {
                        transform: "translate(200, " + f + ")"
                    }));
                    f += 80
                } else {
                    if (!s.subtiers) continue;
                    for (var Q = (s.glyphCacheOrigin - this.viewStart) * this.scale, y = !1, H = 0; H < s.subtiers.length; ++H) {
                        for (var f = f + 3, F = s.subtiers[H], c = [], w = 0; w < F.glyphs.length; ++w) c.push(F.glyphs[w].toSVG());
                        u.appendChild(n(v, "g", c, {
                            transform: "translate(" + (200 +
                                Q) + ", " + f + ")"
                        }));
                        if (F.quant) {
                            var y = !0,
                                c = F.quant,
                                B = F.height,
                                w = 2;
                            40 < B && (w = 1 + (B / 20 | 0));
                            var B = B / (w - 1),
                                I = (c.max - c.min) / (w - 1),
                                x = new z;
                            x.moveTo(205, f);
                            x.lineTo(200, f);
                            x.lineTo(200, f + F.height);
                            x.lineTo(205, f + F.height);
                            for (var C = 1; C < w - 1; ++C) {
                                var A = C * B;
                                x.moveTo(200, f + A);
                                x.lineTo(203, f + A)
                            }
                            G.appendChild(n(v, "path", null, {
                                d: x.toPathData(),
                                fill: "none",
                                stroke: "black",
                                strokeWidth: "2px"
                            }));
                            G.appendChild(n(v, "text", a(c.max), {
                                x: 197,
                                y: f + 7,
                                textAnchor: "end"
                            }));
                            G.appendChild(n(v, "text", a(c.min), {
                                x: 197,
                                y: f + F.height,
                                textAnchor: "end"
                            }));
                            for (C = 1; C < w - 1; ++C) A = C * B, G.appendChild(n(v, "text", a(1 * c.max - C * I), {
                                x: 197,
                                y: f + A + 3,
                                textAnchor: "end"
                            }))
                        }
                        f += F.height + 3
                    }
                    f - L < e && (f = L + e)
                }
                G.appendChild(n(v, "text", "string" === typeof s.config.name ? s.config.name : s.dasSource.name, {
                    x: 200 - (y ? 20 : 12),
                    y: (f + L + 8) / 2,
                    fontSize: "10pt",
                    textAnchor: "end"
                }));
                N.setAttribute("height", f - L);
                d.appendChild(n(v, "g", [u, G]))
            }
            if (b.highlights)
                for (e = this.highlights || [], y = 0; y < e.length; ++y) B = e[y], (B.chr == this.chr || B.chr == "chr" + this.chr) && B.min < this.viewEnd && B.max > this.viewStart && (p = (Math.max(B.min,
                    this.viewStart) - this.viewStart) * this.scale, s = (Math.min(B.max, this.viewEnd) - this.viewStart) * this.scale, d.appendChild(n(v, "rect", null, {
                    x: 200 + p,
                    y: 70,
                    width: s - p,
                    height: f - 70,
                    stroke: "none",
                    fill: this.defaultHighlightFill,
                    fillOpacity: this.defaultHighlightAlpha
                })));
            e = -1;
            "center" == b.ruler ? e = 200 + (this.viewEnd - this.viewStart) * this.scale / 2 : "left" == b.ruler ? e = 200 : "right" == b.ruler && (e = 200 + (this.viewEnd - this.viewStart) * this.scale);
            0 <= e && d.appendChild(n(v, "line", null, {
                x1: e,
                y1: 70,
                x2: e,
                y2: f,
                stroke: "blue"
            }));
            l.appendChild(d);
            g.documentElement.setAttribute("width", this.featurePanelWidth + 20 + 200);
            g.documentElement.setAttribute("height", f + 50);
            return new Blob([(new XMLSerializer).serializeToString(g)], {
                type: "image/svg+xml"
            })
        }
    }, {
        "./cbrowser": 6,
        "./numformats": 26,
        "./sequence-draw": 31,
        "./svg-utils": 38,
        "./utils": 48,
        "./version": 50
    }],
    38: [function(e, u, s) {
        function p() {
            this.ops = []
        }
        p.prototype.moveTo = function(e, m) {
            this.ops.push("M " + e + " " + m)
        };
        p.prototype.lineTo = function(e, m) {
            this.ops.push("L " + e + " " + m)
        };
        p.prototype.closePath = function() {
            this.ops.push("Z")
        };
        p.prototype.toPathData = function() {
            return this.ops.join(" ")
        };
        "undefined" !== typeof u && (u.exports = {
            NS_SVG: "http://www.w3.org/2000/svg",
            NS_XLINK: "http://www.w3.org/1999/xlink",
            SVGPath: p
        })
    }, {}],
    39: [function(e, u, s) {
        function p(b) {
            h.call(this);
            this.readiness = "Connecting";
            this.source = b;
            this.tabixHolder = new z;
            var e = this,
                g = m(b.payload);
            if (g) this.parser = g;
            else throw "Unsuported tabix payload " + b.payload;
            var l;
            this.source.blob ? (b = new r(this.source.blob), l = new r(this.source.indexBlob)) : (b = new v(this.source.uri, {
                credentials: this.source.credentials
            }), l = new v(this.source.indexURI || this.source.uri + ".tbi", {
                credentials: this.source.credentials
            }));
            a(b, l, function(a, b) {
                e.tabixHolder.provide(a);
                a.fetchHeader(function(a, b) {
                    if (a) {
                        for (var d = g.createSession(function() {}), f = 0; f < a.length; ++f) d.parse(a[f]);
                        d.flush()
                    }
                });
                e.readiness = null;
                e.notifyReadiness()
            })
        }
        if ("undefined" !== typeof e) {
            u = e("./sourceadapters");
            var n = u.registerSourceAdapterFactory,
                m = u.makeParser,
                h = u.FeatureSourceBase;
            u = e("./bin");
            var v = u.URLFetchable,
                r = u.BlobFetchable,
                z = e("./utils").Awaited,
                a = e("./tabix").connectTabix
        }
        p.prototype = Object.create(h.prototype);
        p.prototype.fetch = function(a, e, g, h, f, d, m) {
            var n = this;
            n.busy++;
            n.notifyActivity();
            this.tabixHolder.await(function(d) {
                d.fetch(a, e, g, function(a, b) {
                    n.busy--;
                    n.notifyActivity();
                    for (var d = [], f = n.parser.createSession(function(a) {
                            d.push(a)
                        }), g = 0; g < a.length; ++g) f.parse(a[g]);
                    f.flush();
                    m(null, d, 1E9)
                })
            })
        };
        p.prototype.getStyleSheet = function(a) {
            this.parser && this.parser.getStyleSheet && this.parser.getStyleSheet(a)
        };
        p.prototype.getDefaultFIPs =
            function(a) {
                this.parser && this.parser.getDefaultFIPs && this.parser.getDefaultFIPs(a)
            };
        n("tabix", function(a) {
            return {
                features: new p(a)
            }
        })
    }, {
        "./bin": 4,
        "./sourceadapters": 34,
        "./tabix": 40,
        "./utils": 48
    }],
    40: [function(e, u, s) {
        function p() {}

        function n(a, e, g) {
            var l = new p;
            l.data = a;
            l.tbi = e;
            l.tbi.fetch(function(a) {
                if (!a) return g(null, "Couldn't access Tabix");
                a = r(a, a.byteLength);
                var b = new Uint8Array(a);
                if (h(b, 0) != m) return g(null, "Not a tabix index");
                var e = h(b, 4);
                l.format = h(b, 8);
                l.colSeq = h(b, 12);
                l.colStart = h(b, 16);
                l.colEnd = h(b, 20);
                l.meta = h(b, 24);
                l.skip = h(b, 28);
                h(b, 32);
                l.indices = [];
                var n = 36;
                l.chrToIndex = {};
                l.indexToChr = [];
                for (var p = 0; p < e; ++p) {
                    for (var q = "";;) {
                        var s = b[n++];
                        if (0 == s) break;
                        q += String.fromCharCode(s)
                    }
                    l.chrToIndex[q] = p;
                    0 == q.indexOf("chr") ? l.chrToIndex[q.substring(3)] = p : l.chrToIndex["chr" + q] = p;
                    l.indexToChr.push(q)
                }
                q = 1E9;
                for (s = 0; s < e; ++s) {
                    for (var u = n, z = h(b, n), n = n + 4, p = 0; p < z; ++p) {
                        h(b, n);
                        var y = h(b, n + 4),
                            n = n + (8 + 16 * y)
                    }
                    for (var y = h(b, n), H = n += 4, p = 0; p < y; ++p) {
                        var F = v(b, H),
                            H = H + 8;
                        if (F) {
                            p = F.block;
                            0 < F.offset && (p += 65536);
                            p < q && (q = p);
                            break
                        }
                    }
                    n += 8 * y;
                    0 < z && (l.indices[s] = new Uint8Array(a, u, n - u))
                }
                l.headerMax = q;
                g(l)
            })
        }
        var m = 21578324;
        if ("undefined" !== typeof e) {
            e("./spans");
            var h = e("./bin").readInt;
            e = e("./lh3utils");
            var v = e.readVob,
                r = e.unbgzf,
                z = e.reg2bins,
                a = e.Chunk
        }
        p.prototype.blocksForRange = function(b, e, g) {
            var l = this.indices[b];
            if (!l) return [];
            b = z(e, g);
            for (var f = [], d = 0; d < b.length; ++d) f[b[d]] = !0;
            b = [];
            for (var m = [], d = h(l, 0), n = 4, p = 0; p < d; ++p) {
                var r = h(l, n),
                    s = h(l, n + 4),
                    n = n + 8;
                if (f[r])
                    for (var u = 0; u < s; ++u) {
                        var Q = v(l, n),
                            y = v(l, n + 8);
                        (4681 >
                            r ? m : b).push(new a(Q, y));
                        n += 16
                    } else n += 16 * s
            }
            d = h(l, n);
            f = null;
            e = Math.min(e >> 14, d - 1);
            g = Math.min(g >> 14, d - 1);
            for (d = e; d <= g; ++d)(e = v(l, n + 4 + 8 * d)) && (!f || e.block < f.block || e.offset < f.offset) && (f = e);
            l = [];
            if (null != f)
                for (d = 0; d < m.length; ++d) g = m[d], g.maxv.block >= f.block && g.maxv.offset >= f.offset && l.push(g);
            m = l;
            l = [];
            for (d = 0; d < m.length; ++d) l.push(m[d]);
            for (d = 0; d < b.length; ++d) l.push(b[d]);
            l.sort(function(a, b) {
                var c = a.minv.block - b.minv.block;
                return 0 != c ? c : a.minv.offset - b.minv.offset
            });
            b = [];
            if (0 < l.length) {
                m = l[0];
                for (d =
                    1; d < l.length; ++d) g = l[d], g.minv.block == m.maxv.block ? m = new a(m.minv, g.maxv) : (b.push(m), m = g);
                b.push(m)
            }
            return b
        };
        p.prototype.fetch = function(a, e, g, h) {
            function f() {
                if (s >= n.length) return h(p);
                if (v) {
                    var a = new Uint8Array(v);
                    d.readRecords(a, n[s].minv.offset, p, e, g, m);
                    v = null;
                    ++s;
                    return f()
                }
                var b = n[s],
                    a = b.minv.block;
                d.data.slice(a, b.maxv.block + 65536 - a).fetch(function(a) {
                    v = r(a, b.maxv.block - b.minv.block + 1);
                    return f()
                })
            }
            var d = this;
            a = this.chrToIndex[a];
            if (void 0 == a) return h([]);
            var m = this.indexToChr[a],
                n;
            void 0 ===
                a ? n = [] : (n = this.blocksForRange(a, e, g)) || h(null, "Error in index fetch");
            var p = [],
                s = 0,
                v;
            f()
        };
        p.prototype.readRecords = function(a, e, g, h, f, d) {
            a: for (;;) {
                for (var m = ""; e < a.length;) {
                    var n = a[e++];
                    if (10 == n) {
                        n = m.split("\t");
                        if (n[this.colSeq - 1] == d) {
                            var p = parseInt(n[this.colStart - 1]),
                                r = p;
                            0 < this.colEnd && (r = parseInt(n[this.colEnd - 1]));
                            this.format & 65536 && ++p;
                            p <= f && r >= h && g.push(m)
                        }
                        continue a
                    } else m += String.fromCharCode(n)
                }
                break
            }
        };
        p.prototype.fetchHeader = function(a) {
            var e = this;
            e.data.slice(0, e.headerMax).fetch(function(g) {
                if (!g) return a(null,
                    "Fetch failed");
                g = new Uint8Array(r(g, g.byteLength));
                for (var h = 0, f = "", d = []; h < g.length;) {
                    var n = g[h++];
                    if (10 == n)
                        if (f.charCodeAt(0) == e.meta) d.push(f), f = "";
                        else return a(d);
                    else f += String.fromCharCode(n)
                }
                a(d)
            })
        };
        "undefined" !== typeof u && (u.exports = {
            connectTabix: n,
            TABIX_MAGIC: m
        })
    }, {
        "./bin": 4,
        "./lh3utils": 24,
        "./spans": 35
    }],
    41: [function(e, u, s) {
        function p(a) {
            this.genomes = {};
            this.url = a
        }

        function n() {}

        function m(a) {
            this.hub = a
        }

        function h(a, b, f) {
            f = f || {};
            f.salt = !0;
            r(a, function(e, h) {
                    if (h) return b(null, h);
                    var n = e.split(l),
                        q = new p(a);
                    f.credentials && (q.credentials = f.credentials);
                    for (var s = 0; s < n.length - 2; s += 3) q[n[s + 1]] = n[s + 2];
                    if (q.genomesFile) {
                        var y = z(a, q.genomesFile);
                        r(y, function(a, d) {
                            if (d) return b(null, d);
                            for (var c = a.split(g), e = 0; e < c.length; ++e) {
                                var h = c[e].split(l),
                                    n = new m(q);
                                f.credentials && (n.credentials = f.credentials);
                                for (var p = 0; p < h.length - 2; p += 3) n[h[p + 1]] = h[p + 2];
                                n.twoBitPath && (n.twoBitPath = z(y, n.twoBitPath));
                                n.genome && n.trackDb && (n.absURL = z(y, n.trackDb), q.genomes[n.genome] = n)
                            }
                            b(q)
                        }, f)
                    } else b(null, "No genomesFile")
                },
                f)
        }

        function v(a, b) {
            return a.priority && b.priority ? 1 * a.priority - 1 * b.priority : a.priority ? 1 : b.priority ? -1 : a.shortLabel.localeCompare(b.shortLabel)
        }
        if ("undefined" !== typeof e) {
            s = e("./utils");
            var r = s.textXHR,
                z = s.relativeURL,
                a = s.shallowCopy;
            e = e("./das");
            var b = e.DASStylesheet,
                q = e.DASStyle
        }
        var g = /\n\s*\n/,
            l = /(\w+) +(.+)\n?/,
            f = /subGroup[1-9]/;
        n.prototype.get = function(a) {
            if (this[a]) return this[a];
            if (this._parent) return this._parent.get(a)
        };
        m.prototype.getTracks = function(a) {
            var b = this;
            if (this._tracks) return a(this._tracks);
            r(this.absURL, function(e, h) {
                if (h) return a(null, h);
                e = e.replace("\\\n", " ");
                for (var m = [], p = {}, q = e.split(g), r = 0; r < q.length; ++r) {
                    var s = q[r].replace(/\#.*/g, "").split(l),
                        v = new n;
                    v._db = b;
                    for (var u = 0; u < s.length - 2; u += 3) {
                        var c = s[u + 1],
                            w = s[u + 2];
                        if (c.match(f)) {
                            v.subgroups || (v.subgroups = {});
                            for (var c = w.split(/\s/), w = c[0], B = {
                                    name: c[1],
                                    tags: [],
                                    titles: []
                                }, I = 2; I < c.length; ++I) {
                                var x = c[I].split(/=/);
                                B.tags.push(x[0]);
                                B.titles.push(x[1])
                            }
                            v.subgroups[w] = B
                        } else if ("subGroups" === c)
                            for (c = w.split(/(\w+)=(\w+)/), v.sgm = {},
                                I = 0; I < c.length - 2; I += 3) v.sgm[c[I + 1]] = c[I + 2];
                        else v[s[u + 1]] = s[u + 2]
                    }
                    v.track && (v.type || v.container || v.view || v.bigDataUrl) && (m.push(v), p[v.track] = v)
                }
                q = [];
                r = [];
                for (c = 0; c < m.length; ++c) v = m[c], s = !0, v.parent && (u = v.parent.split(/\s+/), (w = p[u[0]]) ? (v._parent = w, w.children || (w.children = []), w.children.push(v), w && (s = !1)) : console.log("Couldn't find parent " + u[0] + "(" + v.parent + ")")), v.compositeTrack ? r.push(v) : s && q.push(v);
                for (m = 0; m < r.length; ++m)
                    if (p = r[m], p.children) {
                        v = !1;
                        for (s = 0; s < p.children.length; ++s) c = p.children[s],
                            c.view && (c.shortLabel = p.shortLabel + ": " + c.shortLabel, q.push(c), v = !0);
                        v || q.push(p)
                    }
                b._tracks = q;
                return a(b._tracks, null)
            }, {
                credentials: this.credentials,
                salt: !0
            })
        };
        n.prototype.toDallianceSource = function() {
            var a = {
                name: this.shortLabel,
                desc: this.longLabel
            };
            this._db.mapping && (a.mapping = this._db.mapping);
            var b = this.get("pennantIcon");
            b && (b = b.split(/\s+/), a.pennant = "http://genome.ucsc.edu/images/" + b[0]);
            if (b = this.get("searchTrix")) a.trixURI = z(this._db.absURL, b);
            if ("multiWig" == this.container) {
                a.merge = "concat";
                a.overlay = [];
                b = this.children || [];
                a.style = [];
                a.noDownsample = !0;
                for (var f = 0; f < b.length; ++f) {
                    var g = b[f],
                        e = g.toDallianceSource();
                    a.overlay.push(e);
                    if (e.style)
                        for (var h = 0; h < e.style.length; ++h) {
                            var l = e.style[h];
                            l.method = g.shortLabel;
                            "transparentOverlay" == this.aggregate && (l.style.ALPHA = 0.5);
                            a.style.push(l)
                        }
                }
                return a
            }
            b = this.type;
            if (!b) {
                for (b = this; b._parent && !b.type;) b = b._parent;
                b = b.type
            }
            if (b) {
                f = b.split(/\s+/);
                if ("bigBed" == f[0] && this.bigDataUrl) return b = f[1] | 0, f = "+" == f[2], a.bwgURI = z(this._db.absURL, this.bigDataUrl),
                    a.style = this.bigbedStyles(), this._db.credentials && (a.credentials = !0), 12 <= b && f && (a.collapseSuperGroups = !0), a;
                if ("bigWig" == f[0] && this.bigDataUrl) return a.bwgURI = z(this._db.absURL, this.bigDataUrl), a.style = this.bigwigStyles(), a.noDownsample = !0, this.yLineOnOff && "on" == this.yLineOnOff && (a.quantLeapThreshold = void 0 !== this.yLineMark ? 1 * this.yLineMark : 0), this._db.credentials && (a.credentials = !0), a;
                if ("bam" == f[0] && this.bigDataUrl) return a.bamURI = z(this._db.absURL, this.bigDataUrl), this._db.credentials && (a.credentials = !0), a;
                if ("vcfTabix" == f[0] && this.bigDataUrl) return a.uri = z(this._db.absURL, this.bigDataUrl), a.tier_type = "tabix", a.payload = "vcf", this._db.credentials && (a.credentials = !0), a;
                console.log("Unsupported " + this.type)
            }
        };
        n.prototype.bigwigStyles = function() {
            var a = this.type;
            if (!a) {
                for (a = this; a._parent && !a.type;) a = a._parent;
                a = a.type
            }
            if (a) {
                var a = a.split(/\s+/),
                    f, g;
                3 <= a.length && (f = 1 * a[1], g = 1 * a[2]);
                var e;
                this.maxHeightPixels && (a = this.maxHeightPixels.split(/:/), 3 == a.length ? e = a[1] | 0 : console.log("maxHeightPixels should be of the form max:default:min"));
                a = "bars";
                this.graphTypeDefault && (a = this.graphTypeDefault);
                var h = "black",
                    l = null;
                this.color && (h = "rgb(" + this.color + ")");
                this.altColor && (l = "rgb(" + this.altColor + ")");
                var n = new b,
                    m = new q;
                m.glyph = "points" == a ? "POINT" : "HISTOGRAM";
                l ? (m.COLOR1 = h, m.COLOR2 = l) : m.BGCOLOR = h;
                m.HEIGHT = e || 30;
                if (f || g) m.MIN = f, m.MAX = g;
                n.pushStyle({
                    type: "default"
                }, null, m);
                return n.styles
            }
        };
        n.prototype.bigbedStyles = function() {
            var d = "on" == ("" + this.get("itemRgb")).toLowerCase(),
                f = this.get("visibility") || "full",
                g = this.get("color"),
                g = g ? "rgb(" +
                g + ")" : "blue",
                e = new b,
                h = new q;
            h.glyph = "BOX";
            h.FGCOLOR = "black";
            h.BGCOLOR = g;
            h.HEIGHT = "full" == f || "pack" == f ? 12 : 8;
            h.BUMP = "full" == f || "pack" == f;
            h.LABEL = "full" == f || "pack" == f;
            h.ZINDEX = 20;
            d && (h.BGITEM = !0);
            (f = this.get("colorByStrand")) ? (f = f.split(/\s+/), g = a(h), g.BGCOLOR = "rgb(" + f[0] + ")", e.pushStyle({
                type: "bigwig",
                orientation: "+"
            }, null, g), h = a(h), h.BGCOLOR = "rgb(" + f[1] + ")", e.pushStyle({
                type: "bigwig",
                orientation: "-"
            }, null, h)) : e.pushStyle({
                type: "bigwig"
            }, null, h);
            h = new q;
            h.glyph = "BOX";
            h.FGCOLOR = "black";
            d && (h.BGITEM = !0);
            h.BGCOLOR = "red";
            h.HEIGHT = 10;
            h.BUMP = !0;
            h.ZINDEX = 20;
            e.pushStyle({
                type: "translation"
            }, null, h);
            d = new q;
            d.glyph = "BOX";
            d.FGCOLOR = "black";
            d.BGCOLOR = "white";
            d.HEIGHT = 10;
            d.ZINDEX = 10;
            d.BUMP = !0;
            d.LABEL = !0;
            e.pushStyle({
                type: "transcript"
            }, null, d);
            return e.styles
        };
        "undefined" !== typeof u && (u.exports = {
            connectTrackHub: h,
            THUB_COMPARE: v
        })
    }, {
        "./das": 10,
        "./utils": 48
    }],
    42: [function(e, u, s) {
        if ("undefined" !== typeof e) var p = e("./cbrowser").Browser,
            n = e("./utils").shallowCopy;
        p.prototype.mergeSelectedTiers = function() {
            for (var e = [], h = [], p = 0; p < this.selectedTiers.length; ++p) {
                var r = this.tiers[this.selectedTiers[p]];
                e.push(n(r.dasSource));
                for (var s = r.stylesheet.styles, a = 0; a < s.length; ++a) {
                    var b = s[a],
                        q = n(b);
                    q.method = r.dasSource.name.replace(/[()+*?]/g, "\\$&");
                    q._methodRE = null;
                    q.style = n(b.style);
                    void 0 === q.style.ZINDEX && (q.style.ZINDEX = p);
                    r.forceMin && (q.style.MIN = r.forceMin);
                    r.forceMax && (q.style.MAX = r.forceMax);
                    h.push(q)
                }
            }
            this.addTier({
                name: "Merged",
                merge: "concat",
                overlay: e,
                noDownsample: !0,
                style: h
            });
            this.setSelectedTier(this.tiers.length -
                1)
        }
    }, {
        "./cbrowser": 6,
        "./utils": 48
    }],
    43: [function(e, u, s) {
        function p(a) {
            for (var b = 0; b < a.styles.length; ++b) {
                var e = a.styles[b].style;
                if ("__SEQUENCE" === e.glyph) return e
            }
        }
        if ("undefined" !== typeof e) {
            var n = e("./cbrowser").Browser,
                m = e("./utils").makeElement;
            u = e("./das");
            var h = u.isDasBooleanTrue,
                v = u.copyStylesheet,
                r = e("./color").dasColourForName
        }
        var z = {
            DOT: !0,
            EX: !0,
            STAR: !0,
            SQUARE: !0,
            CROSS: !0,
            TRIANGLE: !0,
            PLIMSOLL: !0
        };
        n.prototype.openTierPanel = function(a) {
            if ("tier" === this.uiMode && this.manipulatingTier === a) this.hideToolPanel(),
                this.setUiMode("none");
            else {
                var b = function(a) {
                        a.BGGRAD || (1 == w ? ("LINEPLOT" == a.glyph || z[a.glyph] ? a.FGCOLOR = G.value : a.BGCOLOR = G.value, a.COLOR1 = a.COLOR2 = a.COLOR3 = null) : (a.COLOR1 = G.value, a.COLOR2 = L.value, a.COLOR3 = 2 < w ? N.value : null), a._gradient = null, a._plusColor = Q.value, a._minusColor = y.value)
                    },
                    e = function(c) {
                        for (var b = v(a.stylesheet), d = a.browser.zoomForCurrentScale(), f = 0; f < b.styles.length; ++f) {
                            var g = b.styles[f];
                            g.zoom && g.zoom != d || c(g.style)
                        }
                        return b
                    },
                    g = function(c) {
                        a.mergeStylesheet(e(b))
                    };
                this.manipulatingTier =
                    a;
                var l = m("div", null, {
                    className: "tier-edit"
                });
                if (a.dasSource.mapping) {
                    var f = this.chains[a.dasSource.mapping].coords;
                    l.appendChild(m("div", "Mapped from " + f.auth + f.version, null, {
                        background: "gray",
                        paddingBottom: "5px",
                        marginBottom: "5px",
                        textAlign: "center"
                    }))
                }
                var d = m("div", "Editing styles for current zoom level", null, {
                    background: "gray",
                    paddingBottom: "5px",
                    marginBottom: "5px",
                    textAlign: "center",
                    display: "none"
                });
                l.appendChild(d);
                var n = m("input", null, {
                        type: "text"
                    }),
                    s = m("input", null, {
                        type: "checkbox",
                        disabled: this.disablePinning
                    }),
                    u = m("select");
                u.appendChild(m("option", "Histogram", {
                    value: "HISTOGRAM"
                }));
                u.appendChild(m("option", "Line Plot", {
                    value: "LINEPLOT"
                }));
                u.appendChild(m("option", "Ribbon", {
                    value: "GRADIENT"
                }));
                u.appendChild(m("option", "Scatter", {
                    value: "SCATTER"
                }));
                var G = m("input", null, {
                        type: "text",
                        value: "#dd00dd"
                    }),
                    L = m("input", null, {
                        type: "text",
                        value: "#dd00dd"
                    }),
                    N = m("input", null, {
                        type: "text",
                        value: "#dd00dd"
                    }),
                    Q = m("input", null, {
                        type: "text",
                        value: "#ffa07a"
                    }),
                    y = m("input", null, {
                        type: "text",
                        value: "#87cefa"
                    });
                try {
                    G.type = L.type =
                        N.type = "color", Q.type = y.type = "color"
                } catch (H) {}
                var F = [G, L, N],
                    f = m("i", null, {
                        className: "fa fa-plus-circle"
                    }),
                    c = m("i", null, {
                        className: "fa fa-minus-circle"
                    }),
                    w = 1,
                    B = m("td", F),
                    I = function(a) {
                        w = a;
                        for (var c = 0; c < a; ++c) F[c].style.display = "block";
                        for (c = a; c < F.length; ++c) F[c].style.display = "none"
                    };
                f.addEventListener("click", function(a) {
                    3 > w && (I(w + 1), g(null))
                }, !1);
                c.addEventListener("click", function(a) {
                    1 < w && (I(w - 1), g(null))
                }, !1);
                var x = m("input", null, {
                        type: "text",
                        value: "0.0"
                    }),
                    C = m("input", null, {
                        type: "text",
                        value: "10.0"
                    }),
                    A = m("input", null, {
                        type: "checkbox"
                    }),
                    J = m("input", null, {
                        type: "checkbox"
                    }),
                    k = m("input", null, {
                        type: "checkbox",
                        checked: void 0 !== a.quantLeapThreshold
                    }),
                    E = m("input", null, {
                        type: "text",
                        value: a.quantLeapThreshold,
                        disabled: !k.checked
                    }),
                    S = m("input", null, {
                        type: "text",
                        value: "50"
                    }),
                    K = m("input", null, {
                        type: "checkbox"
                    }),
                    R = m("input", null, {
                        type: "text"
                    }),
                    V = m("input", null, {
                        type: "checkbox"
                    }),
                    Z = null;
                0 < a.stylesheet.styles.length && (Z = a.stylesheet.styles[0].style);
                var ja = function() {
                        n.value = "string" === typeof a.config.name ? a.config.name :
                            a.dasSource.name;
                        s.checked = a.pinned;
                        a.forceHeight ? S.value = "" + a.forceHeight : Z && Z.HEIGHT && (S.value = "" + Z.HEIGHT);
                        "number" == typeof a.quantLeapThreshold ? (k.checked = !0, E.disabled = !1, parseFloat(E.value) != a.quantLeapThreshold && (E.value = a.quantLeapThreshold)) : (k.checked = !1, E.disabled = !0);
                        R.value = "number" == typeof a.subtierMax ? "" + a.subtierMax : "" + (a.dasSource.subtierMax || a.browser.defaultSubtierMax);
                        if (0 < a.stylesheet.styles.length) {
                            for (var c = null, b = !1, f = !1, f = a.browser.zoomForCurrentScale(), g = 0, e = 0; e < a.stylesheet.styles.length; ++e) {
                                var l =
                                    a.stylesheet.styles[e];
                                if (!l.zoom || l.zoom == f)
                                    if (++g, l = a.stylesheet.styles[e].style, c || (c = Z = l), "LINEPLOT" == l.glyph || "HISTOGRAM" == l.glyph || "GRADIENT" == l.glyph || h(l.SCATTER)) b || (c = Z = l), b = !0
                            }
                            if (!c) return;
                            d.style.display = g == a.stylesheet.styles.length ? "none" : "block";
                            f = b && 1 == g;
                            b ? (ma.style.display = "table-row", oa.style.display = "table-row", sa.style.display = "none", ta.style.display = "none") : (ma.style.display = "none", oa.style.display = "none", sa.style.display = "table-row", K.checked = h(Z.BUMP), R.disabled = !h(Z.BUMP),
                                ta.style.display = "table-row", V.checked = h(Z.LABEL));
                            f ? (ea.style.display = "table-row", ra.style.display = "table-row") : (ea.style.display = "none", ra.style.display = "none");
                            g = 1;
                            c.COLOR1 ? (G.value = r(c.COLOR1).toHexString(), c.COLOR2 && (L.value = r(c.COLOR2).toHexString(), c.COLOR3 ? (N.value = r(c.COLOR3).toHexString(), g = 3) : g = 2)) : "LINEPLOT" == c.glyph || "DOT" == c.glyph && c.FGCOLOR ? G.value = r(c.FGCOLOR).toHexString() : c.BGCOLOR && (G.value = r(c.BGCOLOR).toHexString());
                            I(g);
                            c._plusColor && (Q.value = r(c._plusColor).toHexString() || c._plusColor);
                            c._minusColor && (y.value = r(c._minusColor).toHexString() || c._minusColor);
                            h(c.SCATTER) ? u.value = "SCATTER" : u.value = c.glyph;
                            var m, w;
                            void 0 !== c.MIN && (g = parseFloat(c.MIN), isNaN(g) || (m = g));
                            a.forceMinDynamic || void 0 === c.MIN && void 0 === a.forceMin ? (A.checked = !1, x.disabled = !0) : (A.checked = !0, x.disabled = !1);
                            void 0 !== c.MAX && (g = parseFloat(c.MAX), isNaN(g) || (w = g));
                            a.forceMaxDynamic || void 0 === c.MAX && void 0 === a.forceMax ? (J.checked = !1, C.disabled = !0) : (J.checked = !0, C.disabled = !1);
                            void 0 != a.forceMin && (m = a.forceMin);
                            void 0 !=
                                a.forceMax && (w = a.forceMax);
                            "number" == typeof m && m != parseFloat(x.value) && (x.value = m);
                            "number" == typeof w && w != parseFloat(C.value) && (C.value = w);
                            (c = p(a.stylesheet)) ? (ca.style.display = "table-row", T.checked = "mismatch" === c.__SEQCOLOR, ka.style.display = "table-row", X.checked = h(c.__INSERTIONS), la.style.display = "table-row", ga.checked = void 0 === c.__disableQuals || !1 === c.__disableQuals, console.log(c.__disableQuals)) : (ca.style.display = "none", ka.style.display = "none", la.style.display = "none");
                            c && T.checked && !f ? ($.style.display =
                                "table-row", pa.style.display = "table-row") : ($.style.display = "none", pa.style.display = "none")
                        }
                        b && a.browser.sourceAdapterIsCapable(a.featureSource, "quantLeap") ? fa.style.display = "table-row" : fa.style.display = "none"
                    },
                    T = m("input", null, {
                        type: "checkbox"
                    }),
                    ca = m("tr", [m("th", "Highlight mismatches & strands"), m("td", T)]);
                T.addEventListener("change", function(c) {
                    c = v(a.stylesheet);
                    p(c).__SEQCOLOR = T.checked ? "mismatch" : "base";
                    a.mergeStylesheet(c)
                });
                var X = m("input", null, {
                        type: "checkbox"
                    }),
                    ka = m("tr", [m("th", "Show insertions"),
                        m("td", X)
                    ]);
                X.addEventListener("change", function(c) {
                    c = v(a.stylesheet);
                    p(c).__INSERTIONS = X.checked ? "yes" : "no";
                    a.mergeStylesheet(c)
                });
                var ga = m("input", null, {
                        type: "checkbox"
                    }),
                    la = m("tr", [m("th", "Reflect base quality as base color transparency"), m("td", ga)]);
                ga.addEventListener("change", function(c) {
                    c = v(a.stylesheet);
                    var b = p(c);
                    b.__disableQuals = !ga.checked;
                    console.log(b.__disableQuals);
                    a.mergeStylesheet(c)
                });
                var ea = m("tr", [m("th", "Style"), m("td", u)]),
                    ra = m("tr", [m("th", ["Colour(s)", f, c]), B]),
                    $ = m("tr", [m("th",
                        "Plus Strand Color"), m("td", Q)]),
                    pa = m("tr", [m("th", "Minus Strand Color"), m("td", y)]),
                    ma = m("tr", [m("th", "Min value"), m("td", [A, " ", x])]),
                    oa = m("tr", [m("th", "Max value"), m("td", [J, " ", C])]),
                    fa = m("tr", [m("th", "Threshold leap:"), m("td", [k, " ", E])]),
                    sa = m("tr", [m("th", "Bump overlaps"), m("td", [K, " limit: ", R])]),
                    ta = m("tr", [m("th", "Label features"), m("td", V)]),
                    f = m("table", [m("tr", [m("th", "Name", {}, {
                            width: "150px",
                            textAlign: "right"
                        }), n]), m("tr", [m("th", "Pin to top"), s]), m("tr", [m("th", "Height"), m("td", S)]),
                        ea, ra, $, pa, ma, oa, fa, sa, ta, ca, ka, la
                    ]);
                ja();
                l.appendChild(f);
                f = m("button", "Reset track", {
                    className: "btn"
                }, {
                    marginLeft: "auto",
                    marginRight: "auto",
                    display: "block"
                });
                f.addEventListener("click", function(c) {
                    a.setConfig({})
                }, !1);
                l.appendChild(f);
                n.addEventListener("input", function(c) {
                    a.mergeConfig({
                        name: n.value
                    })
                }, !1);
                s.addEventListener("change", function(c) {
                    a.mergeConfig({
                        pinned: s.checked
                    })
                }, !1);
                for (f = 0; f < F.length; ++f) F[f].addEventListener("change", g, !1);
                Q.addEventListener("change", g, !1);
                y.addEventListener("change",
                    g, !1);
                u.addEventListener("change", function(c) {
                    c = e(function(a) {
                        "SCATTER" === u.value ? (a.SCATTER = !0, a.glyph = "DOT", a.SIZE = "3") : (a.glyph = u.value, a.SCATTER = void 0);
                        b(a)
                    });
                    a.mergeStylesheet(c)
                }, !1);
                A.addEventListener("change", function(c) {
                    c = {
                        forceMinDynamic: !A.checked
                    };
                    x.disabled = !A.checked;
                    var b = parseFloat(x.value);
                    A.checked && "number" == typeof b && !isNaN(b) && (c.forceMin = parseFloat(b));
                    a.mergeConfig(c)
                });
                x.addEventListener("input", function(c) {
                    c = parseFloat(x.value);
                    "number" != typeof c || isNaN(c) || a.mergeConfig({
                        forceMin: c
                    })
                }, !1);
                J.addEventListener("change", function(c) {
                    c = {
                        forceMaxDynamic: !J.checked
                    };
                    C.disabled = !J.checked;
                    var b = parseFloat(C.value);
                    J.checked && "number" == typeof b && !isNaN(b) && (c.forceMax = parseFloat(b));
                    a.mergeConfig(c)
                });
                C.addEventListener("input", function(c) {
                    c = parseFloat(C.value);
                    "number" != typeof c || isNaN(c) || a.mergeConfig({
                        forceMax: c
                    })
                }, !1);
                S.addEventListener("input", function(c) {
                    c = parseFloat(S.value);
                    "number" != typeof c || isNaN(c) || a.mergeConfig({
                        height: Math.min(500, c | 0)
                    })
                }, !1);
                var U = function() {
                    E.disabled = !k.checked;
                    if (k.checked) {
                        var c = parseFloat(E.value);
                        "number" != typeof c || isNaN(c) || a.mergeConfig({
                            quantLeapThreshold: parseFloat(E.value)
                        })
                    } else a.mergeConfig({
                        quantLeapThreshold: null
                    })
                };
                k.addEventListener("change", function(a) {
                    U()
                }, !1);
                E.addEventListener("input", function(a) {
                    U()
                }, !1);
                V.addEventListener("change", function(c) {
                    c = e(function(a) {
                        a.LABEL = V.checked ? "yes" : "no"
                    });
                    a.mergeStylesheet(c)
                }, !1);
                K.addEventListener("change", function(c) {
                    c = e(function(a) {
                        a.BUMP = K.checked ? "yes" : "no"
                    });
                    a.mergeStylesheet(c)
                }, !1);
                R.addEventListener("input", function(c) {
                    c = parseInt(R.value);
                    "number" == typeof c && 0 < c && a.mergeConfig({
                        subtierMax: c
                    })
                }, !1);
                this.showToolPanel(l);
                this.setUiMode("tier");
                a.addTierListener(ja);
                var Y = a.browser.scale;
                a.browser.addViewListener(function() {
                    a.browser.scale != Y && (Y = a.browser.scale, ja())
                })
            }
        }
    }, {
        "./cbrowser": 6,
        "./color": 9,
        "./das": 10,
        "./utils": 48
    }],
    44: [function(e, u, s) {
        function p(a, b, g, e) {
            this.config = g || {};
            this.id = "tier" + ++q;
            this.browser = a;
            this.dasSource = m(b);
            this.background = e;
            this.viewport = n("canvas",
                null, {
                    width: "" + ((this.browser.featurePanelWidth | 0) + 2E3),
                    height: "30",
                    className: "viewport_12_5"
                }, {
                    position: "inline-block",
                    margin: "0px",
                    border: "0px"
                });
            this.viewportHolder = n("div", this.viewport, {
                className: "viewport-holder_12_5"
            }, {
                background: e,
                position: "absolute",
                padding: "0px",
                margin: "0px",
                border: "0px",
                left: "-1000px",
                minHeight: "200px"
            });
            this.overlay = n("canvas", null, {
                width: +(this.browser.featurePanelWidth | 0),
                height: "30",
                className: "viewport-overlay"
            });
            this.notifier = n("div", "", {
                className: "notifier"
            });
            this.notifierHolder =
                n("div", this.notifier, {
                    className: "notifier-holder"
                });
            this.quantOverlay = n("canvas", null, {
                width: "50",
                height: "56",
                className: "quant-overlay"
            });
            this.removeButton = n("i", null, {
                className: "fa fa-times"
            });
            this.bumpButton = n("i", null, {
                className: "fa fa-plus-circle"
            });
            this.loaderButton = a.makeLoader(16);
            this.loaderButton.style.display = "none";
            this.infoElement = n("div", this.dasSource.desc, {
                className: "track-label-info"
            });
            this.nameButton = n("div", [], {
                className: "tier-tab"
            });
            //this.nameButton.appendChild(this.removeButton);
            b.pennant ? this.nameButton.appendChild(n("img", null, {
                src: b.pennant,
                width: "16",
                height: "16"
            })) : b.mapping && (g = null, this.browser.chains[b.mapping] && (g = this.browser.chains[b.mapping].coords.version), g && this.nameButton.appendChild(n("span", "" + g, null, {
                fontSize: "8pt",
                background: "black",
                color: "white",
                paddingLeft: "3px",
                paddingRight: "3px",
                paddingTop: "1px",
                paddingBottom: "1px",
                marginLeft: "2px",
                borderRadius: "10px"
            })));
            this.nameElement = n("span", b.name);
            this.nameButton.appendChild(n("span", [this.nameElement, this.infoElement], {
                className: "track-name-holder"
            }));
            //this.nameButton.appendChild(this.bumpButton);
            this.nameButton.appendChild(this.loaderButton);
            this.label = n("span", [this.nameButton], {
                className: "btn-group track-label"
            });
            this.row = n("div", [this.viewportHolder, this.overlay, this.quantOverlay], {
                className: "tier" + (b.className ? " " + b.className : "")
            });
            e || (this.row.style.background = "none");
            a.noDefaultLabels || this.row.appendChild(this.label);
            this.row.appendChild(this.notifierHolder);
            this.layoutHeight = 25;
            this.bumped = !0;
            this.styleIdSeed =
                0;
            b.quantLeapThreshold && (this.quantLeapThreshold = b.quantLeapThreshold);
            this.dasSource.collapseSuperGroups && (this.bumped = !1);
            this.layoutWasDone = !1;
            b.featureInfoPlugin && this.addFeatureInfoPlugin(b.featureInfoPlugin);
            this.initSources();
            var h = this;
            this.featureSource && this.featureSource.getDefaultFIPs && !b.noSourceFeatureInfo && this.featureSource.getDefaultFIPs(function(a) {
                a && h.addFeatureInfoPlugin(a)
            });
            this.featureSource && this.featureSource.addReadinessListener && this.featureSource.addReadinessListener(function(a) {
                h.notify(a, -1)
            });
            this.listeners = [];
            this.featuresLoadedListeners = []
        }
        if ("undefined" !== typeof e) {
            s = e("./utils");
            var n = s.makeElement,
                m = s.shallowCopy,
                h = s.miniJSONify;
            s = e("./das");
            var v = s.DASStylesheet,
                r = s.DASStyle,
                z = e("./sha1").b64_sha1;
            s = e("./style");
            var a = s.StyleFilter,
                b = s.StyleFilterSet
        }
        var q = 0;
        p.prototype.setBackground = function(a) {
            this.background = a;
            this.viewportHolder.style.background = a
        };
        p.prototype.toString = function() {
            return this.id
        };
        p.prototype.addFeatureInfoPlugin = function(a) {
            this.featureInfoPlugins || (this.featureInfoPlugins = []);
            this.featureInfoPlugins.push(a)
        };
        p.prototype.init = function() {
            var a = this;
            a.dasSource.style ? (this.setStylesheet({
                styles: a.dasSource.style
            }), this.browser.refreshTier(this)) : (a.status = "Fetching stylesheet", a.fetchStylesheet(function(b, g) {
                if (g || !b) {
                    a.error = "No stylesheet";
                    b = new v;
                    var e = new r;
                    e.glyph = "BOX";
                    e.BGCOLOR = "blue";
                    e.FGCOLOR = "black";
                    b.pushStyle({
                        type: "default"
                    }, null, e);
                    a.setStylesheet(b)
                } else a.setStylesheet(b), b.geneHint && (a.dasSource.collapseSuperGroups = !0, a.bumped = !1, a.updateLabel()), a._updateFromConfig();
                a.browser.refreshTier(a)
            }))
        };
        p.prototype.setStylesheet = function(a) {
            this.baseStylesheet = m(a);
            for (a = 0; a < this.baseStylesheet.styles.length; ++a) {
                var b = this.baseStylesheet.styles[a] = m(this.baseStylesheet.styles[a]);
                b._methodRE = b._labelRE = b._typeRE = null;
                b.style = m(b.style);
                b.style.id = "style" + ++this.styleIdSeed
            }
            this.baseStylesheetValidity = z(h(this.baseStylesheet));
            this._updateFromConfig()
        };
        p.prototype.getSource = function() {
            return this.featureSource
        };
        p.prototype.getDesiredTypes = function(a) {
            if (a = this.getActiveStyleFilters(a)) return a.typeList()
        };
        p.prototype.getActiveStyleFilters = function(f) {
            f = this.browser.zoomForCurrentScale();
            if (this.stylesheet) {
                for (var d = new b, g = this.stylesheet.styles, e = 0; e < g.length; ++e) {
                    var h = g[e];
                    h.zoom && h.zoom != f || d.add(new a(h.type, h.method, h.label))
                }
                return d
            }
        };
        p.prototype.needsSequence = function(a) {
            return this.sequenceSource && 5 > a || (this.dasSource.bamURI || this.dasSource.bamBlob || this.dasSource.bwgURI || this.dasSource.bwgBlob) && 20 > a ? !0 : !1
        };
        p.prototype.viewFeatures = function(a, b, g, e, h) {
            this.currentFeatures = e;
            this.currentSequence =
                h;
            this.notifyFeaturesLoaded();
            this.knownChr = a;
            this.knownCoverage = b;
            this.status && (this.status = null, this._notifierToStatus());
            this.draw()
        };
        p.prototype.draw = function() {
            var a = this.currentSequence;
            this.sequenceSource ? l(this, a) : g(this);
            this.paint();
            this.originHaxx = 0;
            this.browser.arrangeTiers()
        };
        p.prototype.findNextFeature = function(a, b, g, e, h) {
            if (this.quantLeapThreshold) b = b + (this.browser.viewEnd - this.browser.viewStart + 1) * g / 2 | 0, this.featureSource.quantFindNextFeature(a, b, g, this.quantLeapThreshold, h);
            else {
                if (this.knownCoverage &&
                    b >= this.knownCoverage.min() && b <= this.knownCoverage.max() && this.currentFeatures) {
                    for (var l = null, n = 0; n < this.currentFeatures.length; ++n) {
                        var m = this.currentFeatures[n];
                        if (m.min && m.max && !(m.parents && 0 < m.parents.length))
                            if (0 > g)
                                if (1 == e && m.max >= b && m.min < b) {
                                    if (!l || m.min > l.min || m.min == l.min && m.max < l.max) l = m
                                } else m.max < b && (!l || m.max > l.max || m.max == l.max && m.min < l.min || m.min == l.mmin && l.max >= b) && (l = m);
                        else if (1 == e && m.min <= b && m.max > b) {
                            if (!l || m.max < l.max || m.max == l.max && m.min > l.min) l = m
                        } else m.min > b && (!l || m.min < l.min ||
                            m.min == l.min && m.max > l.max || m.max == l.max && l.min <= b) && (l = m)
                    }
                    if (l) return h(l);
                    b = 0 > g ? this.browser.knownSpace.min : this.browser.knownSpace.max
                }
                this.trySourceFNF(a, b, g, h)
            }
        };
        p.prototype.trySourceFNF = function(a, b, g, e) {
            var h = this;
            this.featureSource.findNextFeature(a, b, g, function(a) {
                a || e(a);
                var b = h.browser.getSequenceSource();
                b || e(a);
                b.getSeqInfo(a.segment, function(b) {
                    b ? e(a) : h.trySourceFNF(a.segment, 0 < g ? 1E10 : 0, g, e)
                })
            })
        };
        p.prototype.updateLabel = function() {
            this.bumpButton.className = this.bumped ? "fa fa-minus-circle" :
                "fa fa-plus-circle";
            this.bumpButton.style.display = this.dasSource.collapseSuperGroups ? "inline-block" : "none"
        };
        p.prototype.updateHeight = function() {
            this.currentHeight = Math.max(Math.max(this.layoutHeight, this.label.clientHeight + 2), this.browser.minTierHeight);
            this.row.style.height = "" + this.currentHeight + "px";
            this.browser.updateHeight()
        };
        p.prototype.drawOverlay = function() {
            var a = this.browser,
                b = a.retina && (1 < window.devicePixelRatio || 1 < a.devicePixelRatio);
            this.overlay.height = this.viewport.height;
            this.overlay.width =
                b ? 2 * a.featurePanelWidth : a.featurePanelWidth;
            var g = this.overlay.getContext("2d");
            b && g.scale(2, 2);
            var e = b = a.viewStart,
                h = a.viewEnd;
            if (this.overlayLabelCanvas) {
                var l = (this.glyphCacheOrigin - this.browser.viewStart) * this.browser.scale;
                g.save();
                g.translate(l, 0);
                var m = -l + 2;
                this.dasSource.tierGroup && (m += 15);
                this.overlayLabelCanvas.draw(g, m, a.featurePanelWidth - l);
                g.restore()
            }
            for (l = 0; l < a.highlights.length; ++l) m = a.highlights[l], (m.chr === a.chr || m.chr === "chr" + a.chr) && m.min < h && m.max > e && (g.globalAlpha = a.defaultHighlightAlpha,
                g.fillStyle = a.defaultHighlightFill, g.fillRect((m.min - b) * a.scale, 0, (m.max - m.min) * a.scale, this.overlay.height));
            this.overlay.style.width = a.featurePanelWidth;
            this.overlay.style.height = this.viewport.style.height;
            this.overlay.style.left = "0px"
        };
        p.prototype.updateStatus = function(a) {
            a ? (this.status = a, this._notifierToStatus()) : this.status && (this.status = null, this._notifierToStatus())
        };
        p.prototype.notify = function(a, b) {
            "number" !== typeof b && (b = 2E3);
            this.notifierFadeTimeout && (clearTimeout(this.notifierFadeTimeout),
                this.notifierFadeTimeout = null);
            if (a) {
                if (this._notifierOn(a), 0 < b) {
                    var g = this;
                    this.notifierFadeTimeout = setTimeout(function() {
                        g._notifierToStatus()
                    }, b)
                }
            } else this._notifierToStatus()
        };
        p.prototype._notifierOn = function(a) {
            this.notifier.textContent = a;
            this.notifier.style.opacity = 0.8
        };
        p.prototype._notifierOff = function() {
            this.notifier.style.opacity = 0
        };
        p.prototype._notifierToStatus = function() {
            this.status ? this._notifierOn(this.status) : this._notifierOff()
        };
        p.prototype.setConfig = function(a) {
            this.config = a || {};
            this._updateFromConfig();
            this.notifyTierListeners()
        };
        p.prototype.mergeStylesheet = function(a) {
            this.mergeConfig({
                stylesheet: a,
                stylesheetValidity: this.baseStylesheetValidity
            })
        };
        p.prototype.mergeConfig = function(a) {
            for (var b in a) this.config[b] = a[b];
            this._updateFromConfig();
            this.notifyTierListeners()
        };
        p.prototype._updateFromConfig = function() {
            var a = !1,
                b = !1;
            this.nameElement.textContent = "string" === typeof this.config.name ? this.config.name : this.dasSource.name;
            var g = this.config.height || this.dasSource.forceHeight;
            g != this.forceHeight && (this.forceHeight = g, a = !0);
            this.forceMinDynamic != this.config.forceMinDynamic && (this.forceMinDynamic = this.config.forceMinDynamic, a = !0);
            g = void 0 != this.config.forceMin ? this.config.forceMin : this.dasSource.forceMin;
            this.forceMin != g && (this.forceMin = g, a = !0);
            this.forceMaxDynamic != this.config.forceMaxDynamic && (this.forceMaxDynamic = this.config.forceMaxDynamic, a = !0);
            g = void 0 != this.config.forceMax ? this.config.forceMax : this.dasSource.forceMax;
            this.forceMax != g && (this.forceMax = g, a = !0);
            g = null;
            void 0 !== this.config.quantLeapThreshold ? g = this.config.quantLeapThreshold : void 0 !== this.dasSource.quantLeapThreshold && (g = this.dasSource.quantLeapThreshold);
            g != this.quantLeapThreshold && (this.quantLeapThreshold = g, a = !0);
            g = null;
            this.config.stylesheetValidity == this.baseStylesheetValidity && (g = this.config.stylesheet);
            g = g || this.baseStylesheet;
            this.stylesheet !== g && (this.stylesheet = g, a = !0);
            g = void 0 !== this.config.pinned ? this.config.pinned : this.dasSource.pinned;
            g !== this.pinned && (this.pinned = g, b = !0);
            g = typeof("number" ===
                this.config.subtierMax) ? this.config.subtierMax : this.dasSource.subtierMax || this.browser.defaultSubtierMax;
            g != this.subtierMax && (this.subtierMax = g, a = !0);
            g = void 0 !== this.config.bumped ? this.config.bumped : void 0 !== this.dasSource.bumped ? this.dasSource.bumped : this.dasSource.collapseSuperGroups ? !1 : !0;
            g !== this.bumped && (this.bumped = g, a = !0, this.updateLabel());
            a && this.scheduleRedraw();
            b && this.browser.reorderTiers()
        };
        p.prototype.scheduleRedraw = function() {
            if (this.currentFeatures) {
                var a = this;
                this.redrawTimeout ||
                    (this.redrawTimeout = setTimeout(function() {
                        a.draw();
                        a.redrawTimeout = null
                    }, 10))
            }
        };
        p.prototype.clearTierListeners = function() {
            this.listeners = []
        };
        p.prototype.addTierListener = function(a) {
            this.listeners.push(a)
        };
        p.prototype.notifyTierListeners = function(a) {
            for (var b = 0; b < this.listeners.length; ++b) try {
                this.listeners[b](a)
            } catch (g) {
                console.log(g)
            }
            this.browser.notifyTier()
        };
        p.prototype.clearFeaturesLoadedListeners = function() {
            this.featuresLoadedListeners = []
        };
        p.prototype.addFeaturesLoadedListener = function(a) {
            this.featuresLoadedListeners.push(a)
        };
        p.prototype.notifyFeaturesLoaded = function() {
            for (var a = 0; a < this.featuresLoadedListeners.length; ++a) try {
                this.featuresLoadedListeners[a].call(this)
            } catch (b) {
                console.log(b)
            }
        };
        if ("undefined" !== typeof u) {
            u.exports = {
                DasTier: p
            };
            var g = e("./feature-draw").drawFeatureTier,
                l = e("./sequence-draw").drawSeqTier
        }
    }, {
        "./das": 10,
        "./feature-draw": 18,
        "./sequence-draw": 31,
        "./sha1": 33,
        "./style": 36,
        "./utils": 48
    }],
    45: [function(e, u, s) {
        function p(a, b) {
            for (var d = 0; d < a.length; ++d) {
                var c = a[d];
                c === b ? c.classList.add("active") : c.classList.remove("active")
            }
        }
        if ("undefined" !== typeof e) {
            u = e("./cbrowser");
            var n = u.Browser,
                m = u.sourcesAreEqual;
            u = e("./utils");
            var h = u.makeElement,
                v = u.removeChildren,
                r = u.Observed;
            u = e("./thub");
            var z = u.THUB_COMPARE,
                a = u.connectTrackHub,
                b = e("./domui").makeTreeTableSection,
                q = e("./probe").probeResource;
            u = e("./bin");
            var g = u.URLFetchable,
                l = u.BlobFetchable,
                f = u.readInt,
                d = e("./lh3utils").unbgzf,
                t = e("./bam").BAI_MAGIC,
                D = e("./tabix").TABIX_MAGIC;
            u = e("./das");
            var P = u.DASSource,
                G = u.DASSegment,
                L = u.DASRegistry,
                N = u.coordsMatch,
                Q = e("./encode").EncodeFetchable
        }
        n.prototype.currentlyActive =
            function(a) {
                for (var b = 0; b < this.tiers.length; ++b)
                    if (m(this.tiers[b].dasSource, a)) return this.tiers[b];
                return !1
            };
        n.prototype.makeButton = function(a, b) {
            var d = h("a", a, {
                href: "#"
            });
            b && this.makeTooltip(d, b);
            return h("li", d)
        };
        n.prototype.showTrackAdder = function(e) {
            function m(a) {
                da.style.display = "none";
                ha.style.display = "none";
                ba.style.display = "none";
                U = !1;
                v(M);
                for (var c = h("div", null, {}, {
                        width: "100%"
                    }), d = [], f = 0; f < a.length; ++f) d.push(a[f]);
                d.sort(function(a, c) {
                    return a.shortLabel.toLowerCase().trim().localeCompare(c.shortLabel.toLowerCase().trim())
                });
                a = [];
                for (var g = [], f = 0; f < d.length; ++f) {
                    var e = d[f];
                    e.children && 0 < e.children.length && "multiWig" != e.container ? a.push(e) : g.push(e)
                }
                0 < g.length && a.push({
                    shortLabel: "Others",
                    priority: -1E8,
                    children: g
                });
                a.sort(z);
                for (var l = [], d = 0; d < a.length; ++d) {
                    e = g = a[d];
                    !e.dimensions && e._parent && e._parent.dimensions && (e = e._parent);
                    f = {};
                    if (e.dimensions)
                        for (var n = e.dimensions.split(/(\w+)=(\w+)/), p = 0; p < n.length - 2; p += 3) f[n[p + 1]] = n[p + 2];
                    if (f.dimX && f.dimY) {
                        for (var w = f.dimX, q = f.dimY, f = e.subgroups[w], n = e.subgroups[q], p = {}, r = 0; r <
                            g.children.length; ++r) {
                            var x = g.children[r],
                                e = x.sgm[w],
                                B = x.sgm[q];
                            p[e] || (p[e] = {});
                            p[e][B] = x
                        }
                        q = h("table", null, {
                            className: "table table-striped table-condensed"
                        }, {
                            tableLayout: "fixed"
                        });
                        e = h("tr");
                        e.appendChild(h("th", null, {}, {
                            width: "150px",
                            height: "100px"
                        }));
                        for (B = 0; B < f.titles.length; ++B) w = h("th", h("div", f.titles[B], {}, {
                            transform: "rotate(-60deg)",
                            transformOrigin: "0% 100%",
                            webkitTransform: "rotate(-60deg) translate(20px,10px)",
                            webkitTransformOrigin: "0% 100%",
                            textAlign: "left"
                        }), {}, {
                            width: "35px",
                            height: "100px",
                            verticalAlign: "bottom"
                        }), e.appendChild(w);
                        q.appendChild(e);
                        for (var s = h("tbody", null, {
                                className: "table table-striped table-condensed"
                            }), A = 0; A < n.titles.length; ++A) {
                            var B = n.tags[A],
                                t = h("tr");
                            t.appendChild(h("th", n.titles[A]), {});
                            for (var I = 0; I < f.titles.length; ++I) {
                                var e = f.tags[I],
                                    C = h("td");
                                if (p[e] && p[e][B]) {
                                    e = p[e][B];
                                    r = e.toDallianceSource();
                                    if (!r) continue;
                                    w = h("tr");
                                    x = h("td");
                                    x.style.textAlign = "center";
                                    var u = h("input");
                                    u.type = "checkbox";
                                    u.dalliance_source = r;
                                    qa && (u.dalliance_mapping = qa);
                                    l.push(u);
                                    C.appendChild(u);
                                    u.addEventListener("change", function(a) {
                                        a.target.checked ? k.addTier(a.target.dalliance_source) : k.removeTier(a.target.dalliance_source)
                                    })
                                }
                                t.appendChild(C)
                            }
                            s.appendChild(t)
                        }
                        q.appendChild(s);
                        c.appendChild(b(g.shortLabel, q, 0 == d))
                    } else {
                        n = h("tbody", null, {
                            className: "table table-striped table-condensed"
                        });
                        p = h("table", n, {
                            className: "table table-striped table-condensed"
                        }, {
                            width: "100%",
                            tableLayout: "fixed"
                        });
                        B = 0;
                        g.children.sort(z);
                        for (f = 0; f < g.children.length; ++f)
                            if (e = g.children[f], r = e.toDallianceSource()) w = h("tr"),
                                x = h("td", null, {}, {
                                    width: "30px"
                                }), x.style.textAlign = "center", u = h("input"), u.type = "checkbox", u.dalliance_source = r, qa && (u.dalliance_mapping = qa), l.push(u), x.appendChild(u), u.addEventListener("change", function(a) {
                                    a.target.checked ? k.addTier(a.target.dalliance_source) : k.removeTier(a.target.dalliance_source)
                                }), w.appendChild(x), q = h("td"), q.appendChild(document.createTextNode(e.shortLabel)), e.longLabel && 0 < e.longLabel.length && k.makeTooltip(q, e.longLabel), w.appendChild(q), n.appendChild(w), ++B;
                        1 < a.length || "Others" !==
                            g.shortLabel ? c.appendChild(b(g.shortLabel, p, 0 == d)) : c.appendChild(p)
                    }
                }
                var y = function() {
                    for (var a = 0; a < l.length; ++a) {
                        var c = l[a],
                            b = k.currentlyActive(c.dalliance_source);
                        b ? (c.checked = !0, c.disabled = null != b.sequenceSource) : c.checked = !1
                    }
                };
                y();
                k.addTierListener(function(a) {
                    y()
                });
                M.appendChild(c)
            }

            function n() {
                p(E, la);
                U = "bin";
                da.style.display = "none";
                ha.style.display = "inline";
                ba.style.display = "none";
                v(M);
                var a = h("div", null, {}, {
                    paddingLeft: "10px",
                    paddingRight: "10px"
                });
                
                //a.appendChild(h("h3", "Add custom URL-based data"));
                //a.appendChild(h("p", ["You can add indexed binary data hosted on an web server that supports CORS (", h("a", "full details", {
                //    href: "http://www.biodalliance.org/bin.html"
                //}), ").  Currently supported formats are bigwig, bigbed, and indexed BAM."]));
                
                a.appendChild(h("h3", "Add a Custom Track"));
                a.appendChild(h("p", ["You can add data hosted on another webserver or on your machine.  Supported formats are bigwig, bigbed, and indexed BAM."]));
                
                a.appendChild(h("br"));
                a.appendChild(document.createTextNode("URL: "));
                $ = h("input", "", {
                    size: 80,
                    value: "http://www.biodalliance.org/datasets/ensGene.bb"
                }, {
                    width: "100%"
                });
                a.appendChild($);
                a.appendChild(h("br"));
                a.appendChild(h("b", "- or -"));
                a.appendChild(h("br"));
                a.appendChild(document.createTextNode("File: "));
                fa = h("input", null, {
                    type: "file",
                    multiple: "multiple"
                });
                a.appendChild(fa);
                a.appendChild(h("p", 'Clicking the "Add" button below will initiate a series of test queries.'));
                M.appendChild(a);
                $.focus()
            }

            function c() {
                p(E, ea);
                da.style.display = "none";
                ha.style.display = "inline";
                ba.style.display = "none";
                U = "hub-connect";
                da.style.visibility = "hidden";
                v(M);
                var a = h("div", null, {}, {
                    paddingLeft: "10px",
                    paddingRight: "10px"
                });
                a.appendChild(h("h3", "Connect to a track hub."));
                a.appendChild(h("p", ['Enter the top-level URL (usually points to a file called "hub.txt") of a UCSC-style track hub']));
                $ = h("input", "", {
                    size: 120,
                    value: "http://www.biodalliance.org/datasets/testhub/hub.txt"
                }, {
                    width: "100%"
                });
                a.appendChild($);
                M.appendChild(a);
                $.focus()
            }

            function w() {
                p(E, ga);
                da.style.display = "none";
                ha.style.display = "inline";
                ba.style.display = "none";
                U = "das";
                v(M);
                var a = h("div", null, {}, {
                    paddingLeft: "10px",
                    paddingRight: "10px"
                });
                a.appendChild(h("h3", "Add custom DAS data"));
                a.appendChild(h("p",
                    "This interface is intended for adding custom or lab-specific data.  Public data can be added more easily via the registry interface."));
                a.appendChild(document.createTextNode("URL: "));
                a.appendChild(h("br"));
                $ = h("input", "", {
                    size: 80,
                    value: "http://www.derkholm.net:8080/das/medipseq_reads/"
                }, {
                    width: "100%"
                });
                a.appendChild($);
                a.appendChild(h("p", 'Clicking the "Add" button below will initiate a series of test queries.  If the source is password-protected, you may be prompted to enter credentials.'));
                M.appendChild(a);
                $.focus()
            }

            function B() {
                if (U)
                    if ("das" === U) {
                        var a = $.value.trim();
                        /^.+:\/\//.exec(a) || (a = "http://" + a);
                        a = new P({
                            name: "temporary",
                            uri: a
                        });
                        va(a)
                    } else if ("bin" === U) {
                    var b = fa.files;
                    if (b && 0 < b.length) {
                        a = na = [];
                        U = "multiple";
                        for (var d = 0; d < b.length; ++d) {
                            var f = b[d];
                            f && a.push({
                                blob: f
                            })
                        }
                        for (d = 0; d < a.length; ++d) aa(a[d]);
                        ia()
                    } else a = $.value.trim(), /^.+:\/\//.exec(a) || (a = "http://" + a), b = {
                        uri: a
                    }, a = a.toLowerCase(), 0 == a.indexOf("https://www.encodeproject.org/") && 0 <= a.indexOf("@@download") && (b.transport = "encode"), W(b)
                } else if ("reset" ===
                    U) w();
                else if ("reset-bin" === U) n();
                else if ("reset-hub" === U) c();
                else if ("prompt-bai" === U)(b = fa.files) && 0 < b.length && b[0] ? (Y.baiBlob = b[0], C(Y)) : (b = Y, da.style.display = "none", ha.style.display = "inline", ba.style.display = "inline", v(M), U = "prompt-bai", M.appendChild(h("h2", "Select an index file")), M.appendChild(h("p", "Dalliance requires a BAM index (.bai) file when displaying BAM data.  These normally accompany BAM files.  For security reasons, web applications like Dalliance can only access local files which you have explicity selected.  Please use the file chooser below to select the appropriate BAI file")),
                    M.appendChild(document.createTextNode("Index file: ")), fa = h("input", null, {
                        type: "file"
                    }), M.appendChild(fa), Y = b);
                else if ("prompt-tbi" === U)(b = fa.files) && 0 < b.length && b[0] ? (Y.indexBlob = b[0], A(Y)) : (b = Y, da.style.display = "none", ha.style.display = "inline", ba.style.display = "inline", v(M), U = "prompt-tbi", M.appendChild(h("h2", "Select an index file")), M.appendChild(h("p", "Dalliance requires a Tabix index (.tbi) file when displaying VCF data.  For security reasons, web applications like Dalliance can only access local files which you have explicity selected.  Please use the file chooser below to select the appropriate BAI file")),
                    M.appendChild(document.createTextNode("Index file: ")), fa = h("input", null, {
                        type: "file"
                    }), M.appendChild(fa), Y = b);
                else if ("finalize" === U || "finalize-bin" === U) Y.name = pa.value, b = ma.value, Y.mapping = "__default__" != b ? b : void 0, oa && (Y.maxbins = oa.checked), 1 < sa.value.length && 1 < ta.value.length && (Y.xUser = sa.value, Y.xPass = ta.value), k.addTier(Y), "finalize-bin" == U ? n() : w();
                else if ("hub-connect" === U) a = $.value.trim(), /^.+:\/\//.exec(a) || (a = "http://" + a), s(a);
                else {
                    if ("multiple" === U) {
                        for (b = 0; b < na.length; ++b)
                            if (a = na[b], !a.hidden &&
                                ("bam" != a.tier_type || a.indexBlob || a.indexUri) && ("tabix" != a.tier_type || a.indexBlob || a.indexUri) && (a = wa(a))) a.noPersist = !0, k.addTier(a);
                        n()
                    }
                } else k.removeAllPopups()
            }

            function s(c, b, d) {
                b = b || {};
                for (var f = 0; f < k.hubObjects.length; ++f) {
                    var g = k.hubObjects[f];
                    if (g.hub.url == c) {
                        for (f = 0; f < E.length; ++f) E[f].hub == g && p(E, E[f]);
                        g.getTracks(function(a, c) {
                            c && console.log(c);
                            m(a)
                        });
                        return
                    }
                }
                a(c, function(a, f) {
                    if (f) {
                        if (!d) return s(c, {
                            credentials: !0
                        }, !0);
                        v(M);
                        M.appendChild(h("h2", "Error connecting to track hub"));
                        M.appendChild(h("p",
                            f));
                        U = "reset-hub"
                    } else {
                        var g = null,
                            e = null,
                            l;
                        for (l in a.genomes) {
                            var n = null,
                                w = !1;
                            if (l == k.coordSystem.ucscName) w = !0;
                            else
                                for (var q in k.chains) l == k.chains[q].coords.ucscName && (n = q, w = !0);
                            w && (w = {
                                url: c,
                                genome: l
                            }, b.credentials && (w.credentials = !0), n && (w.mapping = n, a.genomes[l].mapping = n), k.hubs.push(w), k.hubObjects.push(a.genomes[l]), w = ca(a.genomes[l]), ra.appendChild(w), n && g || (g = a.genomes[l], e = w))
                        }
                        g ? (k.notifyTier(), p(E, e), g.getTracks(function(a, c) {
                            m(a)
                        })) : (v(M), M.appendChild(h("h2", "No data for this genome")),
                            M.appendChild(h("p", "This URL appears to be a valid track-hub, but it doesn't contain any data for the coordinate system of this browser")), M.appendChild(h("p", "coordSystem.ucscName = " + k.coordSystem.ucscName)), U = "reset-hub")
                    }
                }, b)
            }

            function x(a, c) {
                var b = a.uri;
                if (c) {
                    var d = /(.+)\/[^\/]+\/?/.exec(b);
                    d && (b = d[1] + "/sources")
                }(new L(b, {
                    credentials: a.credentials
                })).sources(function(b) {
                    if (!b || 0 == b.length) return c ? O(a) : x(a, !0);
                    var d = null;
                    if (1 == b.length) d = b[0];
                    else
                        for (var f = 0; f < b.length; ++f)
                            if (b[f].uri === a.uri) {
                                d =
                                    b[f];
                                break
                            }
                    f = b = !1;
                    if (d && (a.name = d.name, a.desc = d.desc, a.maxbins = d.maxbins ? !0 : !1, d.capabilities && (a.capabilities = d.capabilities), f = !0, d.coords && 1 == d.coords.length))
                        if (d = d.coords[0], N(d, k.coordSystem)) b = !0;
                        else if (k.chains)
                        for (var g in k.chains) N(d, k.chains[g].coords) && (a.mapping = g, b = !0);
                    return O(a, b, f)
                }, function() {
                    return c ? O(a) : x(a, !0)
                })
            }

            function C(a) {
                (a.baiBlob ? new l(a.baiBlob) : "encode" == a.transport ? new Q(a.bamURI + ".bai") : new g(a.bamURI + ".bai", {
                    credentials: a.credentials
                })).slice(0, 256).fetch(function(c) {
                    var b = !1;
                    c && (c = new Uint8Array(c), b = f(c, 0) == t);
                    return b ? O(a, !1, !1, !0) : u("You have selected a valid BAM file, but a corresponding index (.bai) file was not found.  Please index your BAM (samtools index) and place the BAI file in the same directory")
                })
            }

            function A(a) {
                (a.indexBlob ? new l(a.indexBlob) : new g(a.uri + ".tbi")).slice(0, 65536).fetch(function(c) {
                    var b = !1;
                    if (c) {
                        var g = new Uint8Array(c);
                        if (31 == g[0] || 139 == g[1]) c = d(c), g = new Uint8Array(c), b = f(g, 0) == D
                    }
                    return b ? O(a, !1, !1, !0) : u('You have selected a valid VCF file, but a corresponding index (.tbi) file was not found.  Please index your VCF ("tabix -p vcf -f myfile.vcf.gz") and place the .tbi file in the same directory')
                })
            }

            function u(a) {
                da.style.display = "none";
                ha.style.display = "inline";
                ba.style.display = "inline";
                v(M);
                a = a || "Custom data format not recognized";
                M.appendChild(h("h2", "Error adding custom data"));
                M.appendChild(h("p", a));
                M.appendChild(h("p", "Currently supported formats are bigBed, bigWig, and BAM."));
                U = "reset-bin"
            }
            if ("add" === this.uiMode) this.hideToolPanel(), this.setUiMode("none");
            else {
                var k = this;
                e = h("div", null, {
                    className: "dalliance"
                }, {
                    width: "100%",
                    display: "inline-block",
                    boxSizing: "border-box",
                    MozBoxSizing: "border-box",
                    verticalAlign: "top",
                    paddingRight: "15px"
                });
                var E = [],
                    S, K;
                if (!this.noRegistryTabs) {
                    var R = this.makeButton("Registry", "Browse compatible datasources from the DAS registry");
                    E.push(R);
                    for (var V in this.mappableSources)(function(a) {
                        var c = k.makeButton(k.chains[a].srcTag, "Browse datasources mapped from " + k.chains[a].srcTag);
                        E.push(c);
                        c.addEventListener("click", function(b) {
                            b.preventDefault();
                            b.stopPropagation();
                            p(E, c);
                            S(k.mappableSources[a], a)
                        }, !1)
                    })(V)
                }
                V = {};
                for (var Z = 0; Z < this.defaultSources.length; ++Z) {
                    var ja =
                        this.defaultSources[Z],
                        T = ja.group || "Tracks";
                        
                        if (ja['pinned']) { continue; }
                        
                    V[T] ? V[T].push(ja) : V[T] = [ja]
                }
                var ca = function(a) {
                        var b = a.hub,
                            d = h("i", null, {
                                className: "fa fa-list-alt"
                            }, {
                                cursor: "context-menu"
                            }),
                            f = b.shortLabel || "Unknown";
                        a.mapping && (f = f + " (" + a.genome + ")");
                        var f = h("span", [f, " ", d]),
                            g = k.makeButton(f, b.longLabel);
                        g.hub = a;
                        E.push(g);
                        g.addEventListener("click", function(c) {
                            c.preventDefault();
                            c.stopPropagation();
                            p(E, g);
                            v(M);
                            c = k.makeLoader(24);
                            c.style.marginLeft = "auto";
                            c.style.marginRight = "auto";
                            c.style.marginTop = "100px";
                            M.appendChild(h("div",
                                c, null, {
                                    textAlign: "center"
                                }));
                            da.style.display = "none";
                            ha.style.display = "none";
                            ba.style.display = "none";
                            a.getTracks(function(a, c) {
                                c && console.log(c);
                                m(a)
                            })
                        }, !1);
                        d.addEventListener("click", function(b) {
                            b.preventDefault();
                            b.stopPropagation();
                            var d = h("li", h("a", "Remove hub")),
                                f = h("li", h("a", "Enable all")),
                                e = h("li", h("a", "Disable all")),
                                l = h("ul", [d, f, e], {
                                    className: "dropdown-menu"
                                }, {
                                    display: "block"
                                }),
                                m = b.clientX;
                            b = b.clientY;
                            m += document.documentElement.scrollLeft || document.body.scrollLeft;
                            b += document.documentElement.scrollTop ||
                                document.body.scrollTop;
                            l.style.position = "absolute";
                            l.style.top = "" + (b + 10) + "px";
                            l.style.left = "" + (m - 30) + "px";
                            k.hPopupHolder.appendChild(l);
                            var n = function(a) {
                                console.log("cc");
                                document.body.removeEventListener("click", n, !0);
                                k.hPopupHolder.removeChild(l)
                            };
                            document.body.addEventListener("click", n, !0);
                            d.addEventListener("click", function(b) {
                                for (b = 0; b < k.hubObjects.length; ++b)
                                    if (k.hubObjects[b].absURL == a.absURL) {
                                        k.hubObjects.splice(b, 1);
                                        break
                                    }
                                for (b = 0; b < k.hubs.length; ++b) {
                                    var d = k.hubs[b];
                                    "string" === typeof d &&
                                        (d = {
                                            url: d
                                        });
                                    if (d.url == a.hub.url && !d.genome || d.genome == a.genome) {
                                        k.hubs.splice(b, 1);
                                        break
                                    }
                                }
                                k.notifyTier();
                                ra.removeChild(g);
                                p(E, ea);
                                c()
                            }, !1);
                            f.addEventListener("click", function(c) {
                                a.getTracks(function(a, c) {
                                    c && console.log(c);
                                    for (var b = 0; b < a.length; ++b) {
                                        var d = a[b].toDallianceSource();
                                        k.currentlyActive(d) || k.addTier(d)
                                    }
                                })
                            }, !1);
                            e.addEventListener("click", function(c) {
                                a.getTracks(function(a, c) {
                                    c && console.log(c);
                                    for (var b = 0; b < a.length; ++b) {
                                        var d = a[b].toDallianceSource();
                                        k.currentlyActive(d) && k.removeTier(d)
                                    }
                                })
                            }, !1)
                        }, !1);
                        return g
                    },
                    X = null,
                    ka = null;
                for (T in V)(function(a, c) {
                    //var b = k.makeButton(a, "Browse the default set of data for this browser");
                    var b = k.makeButton(a, "");
                    b.addEventListener("click", function(a) {
                        a.preventDefault();
                        a.stopPropagation();
                        p(E, b);
                        S(new r(c))
                    }, !1);
                    E.push(b);
                    X || (X = b, ka = c)
                })(T, V[T]);
                var ga = this.makeButton("DAS", "Add arbitrary DAS data");
                //E.push(ga);
                
                //var la = this.makeButton("Binary", "Add data in bigwig or bigbed format");
                //var la = this.makeButton("Custom", "Add data in bigwig or bigbed format");
                
                var la = this.makeButton("Add a Custom Track", "");
                
                E.push(la);
                for (T = 0; T < this.hubObjects.length; ++T) ca(this.hubObjects[T]);
                var ea = this.makeButton("+",
                    "Connect to a new track-hub");
                //E.push(ea);
                var ra = h("ul", E, {
                    className: "nav nav-tabs"
                }, {
                    marginBottom: "0px"
                });
                e.appendChild(ra);
                var $, pa, ma, oa, fa, sa, ta, U = !1,
                    Y = null,
                    T = h("form", null, {}, {
                        display: "inline-block",
                        width: "100%"
                    });
                T.addEventListener("submit", function(a) {
                    a.stopPropagation();
                    a.preventDefault();
                    B();
                    return !1
                }, !0);
                var M = h("div");
                M.style.position = "relative";
                M.style.overflow = "scroll";
                T.appendChild(M);
                var qa, ua;
                S = function(a, c) {
                    da.style.display = "none";
                    ha.style.display = "none";
                    ba.style.display = "none";
                    ua &&
                        ua.removeListener(K);
                    qa = c;
                    ua = a;
                    ua.addListenerAndFire(K)
                };
                K = function(a) {
                    U = !1;
                    var c = [];
                    v(M);
                    if (a) {
                        for (var b = h("tbody", null, {
                                className: "table table-striped table-condensed"
                            }, {
                                width: "100%"
                            }), d = h("table", b, {
                                className: "table table-striped table-condensed"
                            }, {
                                width: "100%",
                                tableLayout: "fixed"
                            }), f = 0, g = [], e = 0; e < a.length; ++e) g.push(a[e]);
                        g.sort(function(a, c) {
                            return a.name.toLowerCase().trim().localeCompare(c.name.toLowerCase().trim())
                        });
                        for (e = 0; e < g.length; ++e) {
                            a = g[e];
                            var l = h("tr"),
                                m = h("td", null, {}, {
                                    width: "30px"
                                });
                            m.style.textAlign = "center";
                            if (!a.props || a.props.cors) {
                                var n = h("input");
                                n.type = "checkbox";
                                n.dalliance_source = a;
                                qa && (n.dalliance_mapping = qa);
                                m.appendChild(n);
                                c.push(n);
                                n.addEventListener("change", function(a) {
                                    a.target.checked ? k.addTier(a.target.dalliance_source) : k.removeTier(a.target.dalliance_source)
                                })
                            } else m.appendChild(document.createTextNode("!")), k.makeTooltip(m, h("span", ["This data source isn't accessible because it doesn't support ", h("a", "CORS", {
                                href: "http://www.w3.org/TR/cors/"
                            }), "."]));
                            l.appendChild(m);
                            m = h("td");
                            m.appendChild(document.createTextNode(a.name));
                            a.desc && 0 < a.desc.length && k.makeTooltip(m, a.desc);
                            l.appendChild(m);
                            b.appendChild(l);
                            ++f
                        }
                        var p = function() {
                            for (var a = 0; a < c.length; ++a) {
                                var b = c[a];
                                k.currentlyActive(b.dalliance_source) ? b.checked = !0 : b.checked = !1
                            }
                        };
                        p();
                        k.addTierListener(function(a) {
                            p()
                        });
                        M.appendChild(d)
                    } else M.appendChild(h("p", "Dalliance was unable to retrieve data source information from the DAS registry, please try again later"))
                };
                R && R.addEventListener("click", function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    p(E, R);
                    S(k.availableSources)
                }, !1);
                la.addEventListener("click", function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    n()
                }, !1);
                ea.addEventListener("click", function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    c()
                }, !1);
                ga.addEventListener("click", function(a) {
                    a.preventDefault();
                    a.stopPropagation();
                    w()
                }, !1);
                var ha = h("button", "Add", {
                    className: "btn btn-primary"
                });
                ha.addEventListener("click", function(a) {
                    a.stopPropagation();
                    a.preventDefault();
                    B()
                }, !1);
                var va = function(a, c) {
                        var b = k.knownSpace;
                        if (b) {
                            var d = Math.max(b.min, (b.min + b.max - 100) / 2) | 0,
                                b = new G(b.chr, d, Math.min(d + 99, b.max));
                            a.features(b, {}, function(b, d) {
                                if (d) c ? (v(M), M.appendChild(h("h2", "Custom data not found")), M.appendChild(h("p", "DAS uri: " + a.uri + " is not answering features requests")), U = "reset") : (a.credentials = !0, va(a, !0));
                                else {
                                    var f = /\/([^/]+)\/?$/.exec(a.uri);
                                    f && (a.name = f[1]);
                                    x(a)
                                }
                            })
                        } else alert("Can't confirm track-addition to an uninit browser.")
                    },
                    wa = function(a) {
                        var c = {
                            name: a.name
                        };
                        a.credentials && (c.credentials = a.credentials);
                        a.mapping && "__default__" != a.mapping && (c.mapping = a.mapping);
                        a.transport && (c.transport = a.transport);
                        if ("bwg" == a.tier_type) return a.blob ? c.bwgBlob = a.blob : a.uri && (c.bwgURI = a.uri), c;
                        if ("bam" == a.tier_type) return a.blob ? (c.bamBlob = a.blob, c.baiBlob = a.indexBlob) : (c.bamURI = a.uri, c.baiURI = a.indexUri), c;
                        if ("tabix" == a.tier_type) return c.tier_type = "tabix", c.payload = a.payload, a.blob ? (c.blob = a.blob, c.indexBlob = a.indexBlob) : (c.uri = a.uri, c.indexUri = a.indexUri), c;
                        if ("memstore" == a.tier_type) return c.tier_type = "memstore",
                            c.payload = a.payload, a.blob ? c.blob = a.blob : c.uri = a.uri, c
                    },
                    W = function(a) {
                        q(a, function(a, c) {
                            if (c) v(M), M.appendChild(h("h2", "Couldn't access custom data")), M.appendChild(h("p", "" + c)), U = "reset-bin";
                            else {
                                var b = wa(a);
                                return "bam" == a.tier_type ? C(b) : "tabix" == a.tier_type ? A(b) : O(b, !1, !1, !0)
                            }
                        })
                    },
                    O = function(a, c, b, d) {
                        da.style.display = "none";
                        ha.style.display = "inline";
                        ba.style.display = "inline";
                        v(M);
                        M.appendChild(h("h2", "Add custom data: step 2"));
                        M.appendChild(document.createTextNode("Label: "));
                        pa = h("input", "", {
                            value: a.name
                        });
                        M.appendChild(pa);
                        sa = h("input", "");
                        ta = h("input", "");
                        M.appendChild(h("br"));
                        M.appendChild(h("br"));
                        M.appendChild(h("h4", "Coordinate system: "));
                        ma = h("select", null);
                        ma.appendChild(h("option", k.nameForCoordSystem(k.coordSystem), {
                            value: "__default__"
                        }));
                        if (k.chains)
                            for (var f in k.chains) ma.appendChild(h("option", k.nameForCoordSystem(k.chains[f].coords), {
                                value: f
                            }));
                        ma.value = a.mapping || "__default__";
                        M.appendChild(ma);
                        c ? M.appendChild(h("p", "(Based on server response, probably doesn't need changing.)")) : (M.appendChild(h("p", [h("b", "Warning: "), "unable to determine the correct value from server responses.  Please check carefully."])), M.appendChild(h("p", "If you don't see the mapping you're looking for, please contact thomas@biodalliance.org")));
                        d || (M.appendChild(document.createTextNode("Quantitative: ")), oa = h("input", null, {
                            type: "checkbox",
                            checked: !0
                        }), "undefined" !== typeof a.maxbins && (oa.checked = a.maxbins), M.appendChild(oa), b ? M.appendChild(h("p", "(Based on server response, probably doesn't need changing.)")) : M.appendChild(h("p", [h("b", "Warning: "), "unable to determine correct value.  If in doubt, leave checked."])));
                        a.bwgBlob && M.appendChild(h("p", [h("b", "Warning: "), "data added from local file.  Due to the browser security model, the track will disappear if you reload Dalliance."]));
                        pa.focus();
                        U = "bin" === U || "prompt-bai" === U || "prompt-tbi" === U ? "finalize-bin" : "finalize";
                        Y = a
                    },
                    na = null,
                    aa = function(a) {
                        q(a, function(a, c) {
                            c && (a.error = c);
                            for (var b = [], d = {}, f = {}, g = 0; g < na.length; ++g) {
                                var e = na[g];
                                "bam" != e.tier_type || e.indexBlob || (d[e.blob.name] =
                                    e);
                                "tabix" != e.tier_type || e.indexBlob || (f[e.blob.name] = e)
                            }
                            for (g = 0; g < na.length; ++g)
                                if (e = na[g], "bai" === e.tier_type) {
                                    var h = /(.+)\.bai$/.exec(e.blob.name);
                                    h && d[h[1]] && (d[h[1]].indexBlob = e.blob, b.push(g))
                                } else "tabix-index" === e.tier_type && (h = /(.+)\.tbi$/.exec(e.blob.name)) && f[h[1]] && (f[h[1]].indexBlob = e.blob, b.push(g));
                            for (d = b.length - 1; 0 <= d; --d) na.splice(b[d], 1);
                            ia()
                        })
                    },
                    ia = function() {
                        v(M);
                        var a = !1,
                            c = h("table", na.filter(function(a) {
                                return !a.hidden
                            }).map(function(c) {
                                h("tr").appendChild(h("td", c.name || c.blob.name));
                                var b;
                                b = c.error ? h("span", "Error", null, {
                                    color: "red"
                                }) : c.tier_type ? c.payload || c.tier_type : k.makeLoader(16);
                                var d, f = "unknown";
                                "bwg" == c.tier_type || "memstore" == c.tier_type ? f = "okay" : "bam" == c.tier_type ? f = c.indexBlob ? "okay" : "needs-index" : "tabix" == c.tier_type && (f = c.indexBlob ? "okay" : "needs-index");
                                if ("okay" == f) {
                                    d = h("select", null, null, {
                                        width: "150px"
                                    });
                                    d.appendChild(h("option", k.nameForCoordSystem(k.coordSystem), {
                                        value: "__default__"
                                    }));
                                    if (k.chains)
                                        for (var g in k.chains) d.appendChild(h("option", k.nameForCoordSystem(k.chains[g].coords), {
                                            value: g
                                        }));
                                    d.value = c.mapping || "__default__";
                                    d.addEventListener("change", function(a) {
                                        c.mapping = d.value;
                                        console.log(c)
                                    }, !1)
                                } else "needs-index" == f && (d = h("span", "Needs index", {}, {
                                    color: "red"
                                }), a = !0);
                                return h("tr", [h("td", c.name || c.blob.name), h("td", b), h("td", d)])
                            }), {
                                className: "table table-striped table-condensed"
                            });
                        M.appendChild(c);
                        if (a) {
                            M.appendChild(h("p", "Some of these files are missing required index (.bai or .tbi) files.  For security reasons, web applications like Dalliance can only access local files which you have explicity selected.  Please use the file chooser below to select the appropriate index file"));
                            M.appendChild(document.createTextNode("Index file(s): "));
                            var b = h("input", null, {
                                type: "file",
                                multiple: "multiple"
                            });
                            M.appendChild(b);
                            b.addEventListener("change", function(a) {
                                console.log("fileset changed");
                                a = b.files || [];
                                for (var c = 0; c < a.length; ++c) {
                                    var d = a[c];
                                    d && (d = {
                                        blob: d,
                                        hidden: !0
                                    }, na.push(d), aa(d))
                                }
                            }, !1)
                        }
                    },
                    ba = h("button", "Cancel", {
                        className: "btn"
                    });
                ba.addEventListener("click", function(a) {
                    a.stopPropagation();
                    a.preventDefault();
                    "finalize-bin" === U ? n() : w()
                }, !1);
                var da = h("button", "Refresh", {
                    className: "btn"
                });
                da.addEventListener("click", function(a) {
                    a.stopPropagation();
                    a.preventDefault();
                    k.queryRegistry(qa)
                }, !1);
                this.makeTooltip(da, "Click to re-fetch data from the DAS registry");
                V = h("div", [ha, " ", ba, " ", da]);
                V.style.margin = "10px";
                T.appendChild(V);
                e.appendChild(T);
                S(k.availableSources);
                this.showToolPanel(e);
                this.setUiMode("add");
                X && (p(E, X), S(new r(ka)))
            }
        }
    }, {
        "./bam": 1,
        "./bin": 4,
        "./cbrowser": 6,
        "./das": 10,
        "./domui": 11,
        "./encode": 12,
        "./lh3utils": 24,
        "./probe": 28,
        "./tabix": 40,
        "./thub": 41,
        "./utils": 48
    }],
    46: [function(e,
        u, s) {
        function p(e, h, p) {
            h.fetchAsText(function(h) {
                if (!h) return p(null, "Couldn't fetch index-index");
                h = h.split(/(.+)([0-9A-F]{10})\n/);
                for (var s = [], a = [], b = 1; b < h.length; b += 3) s.push(h[b]), a.push(parseInt(h[b + 1], 16));
                return p(new n(s, a, e))
            })
        }

        function n(e, h, n) {
            this.keys = e;
            this.offsets = h;
            this.ix = n
        }
        n.prototype.lookup = function(e, h) {
            for (var n, p = (e + "     ").substring(0, 5).toLowerCase(), s = 0; s < this.keys.length; ++s)
                if (0 > p.localeCompare(this.keys[s])) {
                    n = this.ix.slice(this.offsets[s - 1], this.offsets[s] - this.offsets[s -
                        1]);
                    break
                }
            n || (n = this.ix.slice(this.offsets[this.offsets.length - 1]));
            n.fetchAsText(function(a) {
                a = a.split("\n");
                for (var b = 0; b < a.length; ++b)
                    if (0 == a[b].indexOf(e.toLowerCase() + " ")) return h(a[b].split(" "));
                return h(null)
            })
        };
        "undefined" !== typeof u && (u.exports = {
            connectTrix: p
        })
    }, {}],
    47: [function(e, u, s) {
        function p() {}

        function n(a, f) {
            var d = new p;
            d.data = a;
            d.data.slice(0, 12500).fetch(function(a) {
                if (!a) return f(null, "Couldn't access data");
                var g = new Uint8Array(a);
                a = h(g, 0);
                if (a == b) d.readInt = h;
                else if (a == q) d.readInt =
                    v;
                else return f(null, "Not a .2bit file, magic=0x" + a.toString(16));
                a = d.readInt(g, 4);
                if (0 != a) return f(null, "Unsupported version " + a);
                d.seqCount = d.readInt(g, 8);
                d.seqDict = {};
                var e = 16,
                    l = 0,
                    n = 0,
                    p = function() {
                        for (; l < d.seqCount;) {
                            var a = g[e];
                            if (e + a + 6 >= g.length) return d.data.slice(n + e, 12500).fetch(function(a) {
                                n += e;
                                e = 0;
                                g = new Uint8Array(a);
                                p()
                            });
                            ++e;
                            for (var b = "", h = 1; h <= a; ++h) b += String.fromCharCode(g[e++]);
                            a = d.readInt(g, e);
                            e += 4;
                            d.seqDict[b] = new m(d, a);
                            ++l
                        }
                        return f(d)
                    };
                p()
            })
        }

        function m(a, b) {
            this.tbf = a;
            this.offset =
                b
        }
        if ("undefined" !== typeof e) {
            s = e("./bin");
            var h = s.readInt,
                v = s.readIntBE;
            e = e("./spans");
            var r = e.Range,
                z = e.union,
                a = e.intersection
        }
        var b = 440477507,
            q = 1126646042;
        p.prototype.getSeq = function(a) {
            var b = this.seqDict[a];
            b || (b = this.seqDict["chr" + a]);
            return b
        };
        p.prototype.fetch = function(a, b, d, g) {
            var e = this.getSeq(a);
            if (e) {
                if (d <= b) return g("");
                e.fetch(b, d, g)
            } else return g(null, "Couldn't find " + a)
        };
        m.prototype.init = function(a) {
            if (this.seqOffset) return a();
            var b = this;
            b.tbf.data.slice(b.offset, 8).fetch(function(d) {
                if (!d) return a("Fetch failed");
                d = new Uint8Array(d);
                b._length = b.tbf.readInt(d, 0);
                b.nBlockCnt = b.tbf.readInt(d, 4);
                b.tbf.data.slice(b.offset + 8, 8 * b.nBlockCnt + 4).fetch(function(d) {
                    if (!d) return a("Fetch failed");
                    d = new Uint8Array(d);
                    for (var g = null, e = 0; e < b.nBlockCnt; ++e) var h = b.tbf.readInt(d, 4 * e),
                        n = b.tbf.readInt(d, 4 * (e + b.nBlockCnt)),
                        h = new r(h, h + n - 1),
                        g = g ? z(g, h) : h;
                    b.nBlocks = g;
                    b.mBlockCnt = b.tbf.readInt(d, 8 * b.nBlockCnt);
                    b.seqLength = (b._length + 3) / 4 | 0;
                    b.seqOffset = b.offset + 16 + 8 * (b.nBlockCnt + b.mBlockCnt);
                    return a()
                })
            })
        };
        var g = ["T", "C", "A", "G"];
        m.prototype.fetch = function(b, f, d) {
            --b;
            --f;
            var e = this;
            this.init(function(h) {
                if (h) return d(null, h);
                var n = b >> 2;
                h = f + 3 >> 2;
                if (0 > n || h > e.seqLength) return d("Coordinates out of bounds: " + b + ":" + f);
                e.tbf.data.slice(e.seqOffset + n, h - n).salted().fetch(function(h) {
                    function m(a) {
                        for (; u <= a;) {
                            var b = u & 3,
                                d = p[(u >> 2) - n];
                            s += g[0 == b ? d >> 6 & 3 : 1 == b ? d >> 4 & 3 : 2 == b ? d >> 2 & 3 : d & 3];
                            ++u
                        }
                    }
                    if (null == h) return d("SeqFetch failed");
                    var p = new Uint8Array(h);
                    h = [];
                    if (e.nBlocks) {
                        var q = a(new r(b, f), e.nBlocks);
                        q && (h = q.ranges())
                    }
                    for (var s = "", u = b, q = 0; q <
                        h.length; ++q) {
                        var v = h[q];
                        if (u > v.min()) throw "N mismatch...";
                        for (u < v.min() && m(v.min() - 1); u <= v.max();) s += "N", ++u
                    }
                    u <= f && m(f);
                    return d(s)
                })
            })
        };
        m.prototype.length = function(a) {
            var b = this;
            this.init(function(d) {
                return d ? a(null, d) : a(b._length)
            })
        };
        "undefined" !== typeof u && (u.exports = {
            makeTwoBit: n
        })
    }, {
        "./bin": 4,
        "./spans": 35
    }],
    48: [function(e, u, s) {
        function p(a, b) {
            for (var d = 0; d < a.length; ++d)
                if (a[d] == b) return;
            a.push(b)
        }

        function n(a, b, d) {
            a[b] ? a[b].push(d) : a[b] = [d]
        }

        function m(a, b, d) {
            var f = a[b];
            if (f) {
                for (a = 0; a < f.length; ++a)
                    if (f[a] ==
                        d) return;
                f.push(d)
            } else a[b] = [d]
        }

        function h(a, b, d, f) {
            if (a) return a;
            if (b) return b;
            if (d) return d;
            if (f) return f
        }

        function p(a, b) {
            for (var d = 0; d < a.length; ++d)
                if (a[d] == b) return;
            a.push(b)
        }

        function v(a, b) {
            if (!a) return -1;
            for (var d = 0; d < a.length; ++d)
                if (a[d] === b) return d;
            return -1
        }

        function r(a, b, d, f) {
            a = document.createElement(a);
            if (b) {
                b instanceof Array || (b = [b]);
                for (var g = 0; g < b.length; ++g) {
                    var e = b[g];
                    e && ("string" == typeof e ? e = document.createTextNode(e) : "number" == typeof e && (e = document.createTextNode("" + e)), a.appendChild(e))
                }
            }
            if (d)
                for (var c in d) try {
                    a[c] =
                        d[c]
                } catch (h) {
                    throw console.log("error setting " + c), h;
                }
            if (f)
                for (c in f) a.style[c] = f[c];
            return a
        }

        function z(a, b, d, f) {
            a = document.createElementNS(a, b);
            if (d)
                for (d instanceof Array || (d = [d]), b = 0; b < d.length; ++b) {
                    var g = d[b];
                    "string" == typeof g && (g = document.createTextNode(g));
                    a.appendChild(g)
                }
            if (f)
                for (var e in f) {
                    d = a;
                    b = e;
                    var g = f[e],
                        c = P[b];
                    if (!c) {
                        for (var c = "", h = 0; h < b.length; ++h) var l = b.substring(h, h + 1),
                            n = l.toLowerCase(),
                            c = n != l ? c + "-" + n : c + l;
                        P[b] = c
                    }
                    d.setAttribute(c, g)
                }
            return a
        }

        function a(a) {
            if (a && a.childNodes)
                for (; 0 <
                    a.childNodes.length;) a.removeChild(a.firstChild)
        }

        function b(a, d) {
            if ("undefined" === typeof a) return "undefined";
            if (null == a) return "null";
            if ("string" == typeof a) return "'" + a + "'";
            if ("number" == typeof a || "boolean" == typeof a) return "" + a;
            if ("object" == typeof a) {
                if (a instanceof Array) {
                    for (var f = null, g = 0; g < a.length; ++g) f = (null == f ? "" : f + ", ") + b(a[g], d);
                    return "[" + (f ? f : "") + "]"
                }
                d = d || {};
                f = null;
                for (g in a) d[g] || void 0 != g && "function" != typeof a[g] && (f = (null == f ? "" : f + ", ") + g + ": " + b(a[g], d));
                return "{" + (f ? f : "") + "}"
            }
            return typeof a
        }

        function q(a) {
            var b = {},
                d;
            for (d in a) b[d] = a[d];
            return b
        }

        function g(a) {
            this.value = a;
            this.listeners = []
        }

        function l() {
            this.queue = []
        }

        function f(a, b, d) {
            d && d.salt && (a = a + "?salt=" + D("" + Date.now() + "," + ++G));
            var f = new XMLHttpRequest;
            f.onreadystatechange = function() {
                4 == f.readyState && (300 <= f.status ? b(null, "Error code " + f.status) : b(f.responseText))
            };
            f.open("GET", a, !0);
            f.responseType = "text";
            d && d.credentials && (f.withCredentials = !0);
            f.send("")
        }

        function d(a, b) {
            if (0 == b.indexOf("http:") || 0 == b.indexOf("https:")) return b;
            var d = a.lastIndexOf("/");
            return 0 <= d ? a.substr(0, d + 1) + b : b
        }

        function t(a) {
            return r("a", null, {
                href: a
            }).href
        }
        if ("undefined" !== typeof e) var D = e("./sha1").b64_sha1;
        var P = {};
        g.prototype.addListener = function(a) {
            this.listeners.push(a)
        };
        g.prototype.addListenerAndFire = function(a) {
            this.listeners.push(a);
            a(this.value)
        };
        g.prototype.removeListener = function(a) {
            var b = this.listeners;
            a = v(b, a);
            0 <= a && b.splice(a, 1)
        };
        g.prototype.get = function() {
            return this.value
        };
        g.prototype.set = function(a) {
            this.value = a;
            for (var b = 0; b < this.listeners.length; ++b) this.listeners[b](a)
        };
        l.prototype.provide = function(a) {
            if (void 0 !== this.res) throw "Resource has already been provided.";
            this.res = a;
            for (var b = 0; b < this.queue.length; ++b) this.queue[b](a);
            this.queue = null
        };
        l.prototype.await = function(a) {
            if (void 0 !== this.res) return a(this.res), this.res;
            this.queue.push(a)
        };
        var G = 0;
        e = {
            TTT: "F",
            TTC: "F",
            TTA: "L",
            TTG: "L",
            CTT: "L",
            CTC: "L",
            CTA: "L",
            CTG: "L",
            ATT: "I",
            ATC: "I",
            ATA: "I",
            ATG: "M",
            GTT: "V",
            GTC: "V",
            GTA: "V",
            GTG: "V",
            TCT: "S",
            TCC: "S",
            TCA: "S",
            TCG: "S",
            CCT: "P",
            CCC: "P",
            CCA: "P",
            CCG: "P",
            ACT: "T",
            ACC: "T",
            ACA: "T",
            ACG: "T",
            GCT: "A",
            GCC: "A",
            GCA: "A",
            GCG: "A",
            TAT: "Y",
            TAC: "Y",
            TAA: "*",
            TAG: "*",
            CAT: "H",
            CAC: "H",
            CAA: "Q",
            CAG: "Q",
            AAT: "N",
            AAC: "N",
            AAA: "K",
            AAG: "K",
            GAT: "D",
            GAC: "D",
            GAA: "E",
            GAG: "E",
            TGT: "C",
            TGC: "C",
            TGA: "*",
            TGG: "W",
            CGT: "R",
            CGC: "R",
            CGA: "R",
            CGG: "R",
            AGT: "S",
            AGC: "S",
            AGA: "R",
            AGG: "R",
            GGT: "G",
            GGC: "G",
            GGA: "G",
            GGG: "G"
        };
        "trim" in String.prototype || (String.prototype.trim = function() {
            return this.replace(/^\s+/, "").replace(/\s+$/, "")
        });
        "undefined" !== typeof u && (u.exports = {
            textXHR: f,
            relativeURL: d,
            resolveUrlToPage: t,
            shallowCopy: q,
            pusho: n,
            pushnew: p,
            pushnewo: m,
            arrayIndexOf: v,
            pick: h,
            makeElement: r,
            makeElementNS: z,
            removeChildren: a,
            miniJSONify: b,
            Observed: g,
            Awaited: l,
            AMINO_ACID_TRANSLATION: e
        })
    }, {
        "./sha1": 33
    }],
    49: [function(e, u, s) {
        function p() {
            this.info = []
        }

        function n(a, b) {
            this.parser = a;
            this.sink = b
        }
        if ("undefined" !== typeof e) {
            var m = e("./sourceadapters").registerParserFactory;
            e = e("./das");
            var h = e.DASStylesheet,
                v = e.DASStyle,
                r = e.DASFeature
        }
        var z = /([^;=]+)(=([^;]+))?;?/,
            a = /##INFO=<([^>]+)>/,
            b = /([^,=]+)=([^,]+|"[^"]+"),?/;
        p.prototype.createSession =
            function(a) {
                return new n(this, a)
            };
        n.prototype.parse = function(e) {
            if (0 != e.length)
                if ("#" == e[0]) {
                    if (1 < e.length && "#" == e[1] && (e = a.exec(e))) {
                        var g = e[1].split(b),
                            h = null,
                            f = null;
                        for (e = 0; e < g.length - 1; e += 3) {
                            var d = g[e + 1],
                                n = g[e + 2].replace(/"/g, "");
                            "ID" == d ? h = n : "Description" == d && (f = n)
                        }
                        h && f && this.parser.info.push({
                            id: h,
                            desc: f
                        })
                    }
                } else {
                    g = e.split("\t");
                    h = new r;
                    h.segment = g[0];
                    h.id = g[2];
                    h.refAllele = g[3];
                    h.altAlleles = g[4].split(",");
                    h.min = parseInt(g[1]);
                    h.max = h.min + h.refAllele.length - 1;
                    g = g[7].split(z);
                    h.info = {};
                    for (e = 0; e <
                        g.length; e += 4) h.info[g[e + 1]] = g[e + 3];
                    e = h.altAlleles[0];
                    g = h.refAllele;
                    e.length > g.length ? (h.type = "insertion", 0 == e.indexOf(g) ? (h.insertion = e.substr(g.length), h.min += g.length, h.max = h.min - 1) : h.insertion = e) : h.type = e.length < g.length ? "deletion" : "substitution";
                    this.sink(h)
                }
        };
        n.prototype.flush = function() {};
        p.prototype.getStyleSheet = function(a) {
            var b = new h,
                e = new v;
            e.glyph = "__INSERTION";
            e.BUMP = "yes";
            e.LABEL = "no";
            e.FGCOLOR = "rgb(50,80,255)";
            e.BGCOLOR = "#888888";
            e.STROKECOLOR = "black";
            b.pushStyle({
                    type: "insertion"
                },
                null, e);
            e = new v;
            e.glyph = "PLIMSOLL";
            e.BUMP = "yes";
            e.LABEL = "no";
            e.FGCOLOR = "rgb(255, 60, 60)";
            e.BGCOLOR = "#888888";
            e.STROKECOLOR = "black";
            b.pushStyle({
                type: "deletion"
            }, null, e);
            e = new v;
            e.glyph = "PLIMSOLL";
            e.BUMP = "yes";
            e.LABEL = "no";
            e.FGCOLOR = "rgb(50,80,255)";
            e.BGCOLOR = "#888888";
            e.STROKECOLOR = "black";
            b.pushStyle({
                type: "default"
            }, null, e);
            return a(b)
        };
        p.prototype.getDefaultFIPs = function(a) {
            var b = this;
            a(function(a, f) {
                f.add("Ref. allele", a.refAllele);
                f.add("Alt. alleles", a.altAlleles.join(","));
                if (a.info)
                    for (var d =
                            0; d < b.info.length; ++d) {
                        var e = b.info[d],
                            h = a.info[e.id];
                        void 0 !== h && f.add(e.desc, h)
                    }
            })
        };
        m("vcf", function() {
            return new p
        })
    }, {
        "./das": 10,
        "./sourceadapters": 34
    }],
    50: [function(e, u, s) {
        e = {
            CONFIG: 5,
            MAJOR: 0,
            MINOR: 13,
            MICRO: 2,
            PATCH: "",
            BRANCH: "",
            toString: function() {
                var e = "" + this.MAJOR + "." + this.MINOR + "." + this.MICRO;
                this.PATCH && (e += this.PATCH);
                this.BRANCH && "" != this.BRANCH && (e = e + "-" + this.BRANCH);
                return e
            }
        };
        "undefined" !== typeof u && (u.exports = e)
    }, {}],
    51: [function(e, u, s) {
        function p() {
            function e(a) {
                a = Math.min(a, q);
                f =
                    a = Math.max(a, b);
                r.style.left = "" + f + "px"
            }

            function h(a) {
                a = Math.min(a, q);
                d = a = Math.max(a, b);
                s.style.left = "" + d + "px"
            }
            var p = n("hr", null, {
                    className: "slider-track"
                }),
                r = n("hr", null, {
                    className: "slider-thumb active"
                }),
                s = n("hr", null, {
                    className: "slider-thumb"
                }),
                a = n("div", [p, r, s], {
                    className: "slider"
                }),
                b = 0,
                q = 200,
                g = 0,
                l = 200,
                f = 50,
                d = 100,
                u = [];
            a.removeLabels = function() {
                for (var b = 0; b < u.length; ++b) a.removeChild(u[b]);
                u = []
            };
            a.addLabel = function(d, f) {
                var e = n("div", f, {
                    className: "slider-label"
                }, {
                    left: "" + (b + (d - g) * (q - b) / (l - g) | 0) +
                        "px"
                });
                a.appendChild(e);
                u.push(e)
            };
            var D = document.createEvent("HTMLEvents");
            D.initEvent("change", !0, !1);
            Object.defineProperty(a, "value", {
                get: function() {
                    return g + (f - b) * (l - g) / (q - b)
                },
                set: function(a) {
                    e(b + (a - g) * (q - b) / (l - g))
                }
            });
            Object.defineProperty(a, "value2", {
                get: function() {
                    return g + (d - b) * (l - g) / (q - b)
                },
                set: function(a) {
                    h(b + (a - g) * (q - b) / (l - g))
                }
            });
            Object.defineProperty(a, "active", {
                get: function() {
                    return r.classList.contains("active") ? 1 : 2
                },
                set: function(a) {
                    1 == a ? (r.classList.add("active"), s.classList.remove("active")) :
                        (s.classList.add("active"), r.classList.remove("active"))
                }
            });
            Object.defineProperty(a, "min", {
                get: function() {
                    return g
                },
                set: function(a) {
                    g = a
                }
            });
            Object.defineProperty(a, "max", {
                get: function() {
                    return l
                },
                set: function(a) {
                    l = a
                }
            });
            var P, G, p = function(b) {
                G = this == r ? 1 : 2;
                G != a.active && (a.active = G, a.dispatchEvent(D));
                b.stopPropagation();
                b.preventDefault();
                window.addEventListener("mousemove", L, !1);
                window.addEventListener("mouseup", N, !1);
                P = b.clientX - (1 == G ? f : d)
            };
            r.addEventListener("mousedown", p, !1);
            s.addEventListener("mousedown",
                p, !1);
            var L = function(b) {
                    1 == G ? e(b.clientX - P) : h(b.clientX - P);
                    a.dispatchEvent(D)
                },
                N = function(a) {
                    window.removeEventListener("mousemove", L, !1);
                    window.removeEventListener("mouseup", N, !1)
                };
            return a
        }
        if ("undefined" !== typeof e) var n = e("./utils").makeElement;
        "undefined" !== typeof u && (u.exports = p)
    }, {
        "./utils": 48
    }],
    52: [function(e, u, s) {
        u = e("./promise/promise").Promise;
        e = e("./promise/polyfill").polyfill;
        s.Promise = u;
        s.polyfill = e
    }, {
        "./promise/polyfill": 57,
        "./promise/promise": 58
    }],
    53: [function(e, u, s) {
        var p = e("./utils").isArray,
            n = e("./utils").isFunction;
        s.all = function(e) {
            if (!p(e)) throw new TypeError("You must pass an array to all.");
            return new this(function(h, p) {
                function r(b) {
                    return function(e) {
                        s[b] = e;
                        0 === --a && h(s)
                    }
                }
                var s = [],
                    a = e.length,
                    b;
                0 === a && h([]);
                for (var q = 0; q < e.length; q++)(b = e[q]) && n(b.then) ? b.then(r(q), p) : (s[q] = b, 0 === --a && h(s))
            })
        }
    }, {
        "./utils": 62
    }],
    54: [function(e, u, s) {
        (function(e, n) {
            function m() {
                return function() {
                    e.nextTick(r)
                }
            }

            function h() {
                var b = 0,
                    f = new a(r),
                    d = document.createTextNode("");
                f.observe(d, {
                    characterData: !0
                });
                return function() {
                    d.data = b = ++b % 2
                }
            }

            function u() {
                return function() {
                    b.setTimeout(r, 1)
                }
            }

            function r() {
                for (var a = 0; a < q.length; a++) {
                    var b = q[a];
                    (0, b[0])(b[1])
                }
                q = []
            }
            var z = "undefined" !== typeof window ? window : {},
                a = z.MutationObserver || z.WebKitMutationObserver,
                b = "undefined" !== typeof n ? n : void 0 === this ? window : this,
                q = [],
                g;
            g = "undefined" !== typeof e && "[object process]" === {}.toString.call(e) ? m() : a ? h() : u();
            s.asap = function(a, b) {
                1 === q.push([a, b]) && g()
            }
        }).call(this, e("1YiZ5S"), "undefined" !== typeof self ? self : "undefined" !==
            typeof window ? window : {})
    }, {
        "1YiZ5S": 63
    }],
    55: [function(e, u, s) {
        s.cast = function(e) {
            return e && "object" === typeof e && e.constructor === this ? e : new this(function(n) {
                n(e)
            })
        }
    }, {}],
    56: [function(e, u, s) {
        var p = {
            instrument: !1
        };
        s.config = p;
        s.configure = function(e, m) {
            if (2 === arguments.length) p[e] = m;
            else return p[e]
        }
    }, {}],
    57: [function(e, u, s) {
        (function(p) {
            var n = e("./promise").Promise,
                m = e("./utils").isFunction;
            s.polyfill = function() {
                var e;
                e = "undefined" !== typeof p ? p : "undefined" !== typeof window && window.document ? window : self;
                "Promise" in e && "cast" in e.Promise && "resolve" in e.Promise && "reject" in e.Promise && "all" in e.Promise && "race" in e.Promise && function() {
                    var n;
                    new e.Promise(function(e) {
                        n = e
                    });
                    return m(n)
                }() || (e.Promise = n)
            }
        }).call(this, "undefined" !== typeof self ? self : "undefined" !== typeof window ? window : {})
    }, {
        "./promise": 58,
        "./utils": 62
    }],
    58: [function(e, u, s) {
        function p(a) {
            if (!d(a)) throw new TypeError("You must pass a resolver function as the first argument to the promise constructor");
            if (!(this instanceof p)) throw new TypeError("Failed to construct 'Promise': Please use the 'new' operator, this object constructor cannot be called as a function.");
            this._subscribers = [];
            n(a, this)
        }

        function n(a, d) {
            function c(a) {
                z(d, a)
            }

            function f(a) {
                b(d, a)
            }
            try {
                a(c, f)
            } catch (e) {
                f(e)
            }
        }

        function m(a, f, c, e) {
            var g = d(c),
                h, l, n, m;
            if (g) try {
                h = c(e), n = !0
            } catch (p) {
                m = !0, l = p
            } else h = e, n = !0;
            r(f, h) || (g && n ? z(f, h) : m ? b(f, l) : a === Q ? z(f, h) : a === y && b(f, h))
        }

        function h(a, b, c, d) {
            a = a._subscribers;
            var f = a.length;
            a[f] = b;
            a[f + Q] = c;
            a[f + y] = d
        }

        function v(a, b) {
            for (var c, d, f = a._subscribers, e = a._detail, g = 0; g < f.length; g += 3) c = f[g], d = f[g + b], m(b, c, d, e);
            a._subscribers = null
        }

        function r(e, g) {
            var c = null,
                h;
            try {
                if (e ===
                    g) throw new TypeError("A promises callback cannot return that same promise.");
                if (f(g) && (c = g.then, d(c))) return c.call(g, function(c) {
                    if (h) return !0;
                    h = !0;
                    g !== c ? z(e, c) : a(e, c)
                }, function(a) {
                    if (h) return !0;
                    h = !0;
                    b(e, a)
                }), !0
            } catch (l) {
                if (h) return !0;
                b(e, l);
                return !0
            }
            return !1
        }

        function z(b, d) {
            b === d ? a(b, d) : r(b, d) || a(b, d)
        }

        function a(a, b) {
            a._state === L && (a._state = N, a._detail = b, l.async(q, a))
        }

        function b(a, b) {
            a._state === L && (a._state = N, a._detail = b, l.async(g, a))
        }

        function q(a) {
            v(a, a._state = Q)
        }

        function g(a) {
            v(a, a._state = y)
        }
        var l = e("./config").config;
        e("./config");
        var f = e("./utils").objectOrFunction,
            d = e("./utils").isFunction;
        e("./utils");
        u = e("./cast").cast;
        var t = e("./all").all,
            D = e("./race").race,
            P = e("./resolve").resolve,
            G = e("./reject").reject;
        e = e("./asap").asap;
        l.async = e;
        var L = void 0,
            N = 0,
            Q = 1,
            y = 2;
        p.prototype = {
            constructor: p,
            _state: void 0,
            _detail: void 0,
            _subscribers: void 0,
            then: function(a, b) {
                var c = this,
                    d = new this.constructor(function() {});
                if (this._state) {
                    var f = arguments;
                    l.async(function() {
                        m(c._state, d, f[c._state - 1], c._detail)
                    })
                } else h(this,
                    d, a, b);
                return d
            },
            "catch": function(a) {
                return this.then(null, a)
            }
        };
        p.all = t;
        p.cast = u;
        p.race = D;
        p.resolve = P;
        p.reject = G;
        s.Promise = p
    }, {
        "./all": 53,
        "./asap": 54,
        "./cast": 55,
        "./config": 56,
        "./race": 59,
        "./reject": 60,
        "./resolve": 61,
        "./utils": 62
    }],
    59: [function(e, u, s) {
        var p = e("./utils").isArray;
        s.race = function(e) {
            if (!p(e)) throw new TypeError("You must pass an array to race.");
            return new this(function(m, h) {
                for (var p, r = 0; r < e.length; r++)(p = e[r]) && "function" === typeof p.then ? p.then(m, h) : m(p)
            })
        }
    }, {
        "./utils": 62
    }],
    60: [function(e,
        u, s) {
        s.reject = function(e) {
            return new this(function(n, m) {
                m(e)
            })
        }
    }, {}],
    61: [function(e, u, s) {
        s.resolve = function(e) {
            return new this(function(n, m) {
                n(e)
            })
        }
    }, {}],
    62: [function(e, u, s) {
        function p(e) {
            return "function" === typeof e
        }
        e = Date.now || function() {
            return (new Date).getTime()
        };
        s.objectOrFunction = function(e) {
            return p(e) || "object" === typeof e && null !== e
        };
        s.isFunction = p;
        s.isArray = function(e) {
            return "[object Array]" === Object.prototype.toString.call(e)
        };
        s.now = e
    }, {}],
    63: [function(e, u, s) {
        function p() {}
        e = u.exports = {};
        e.nextTick = function() {
            if ("undefined" !== typeof window && window.setImmediate) return function(e) {
                return window.setImmediate(e)
            };
            if ("undefined" !== typeof window && window.postMessage && window.addEventListener) {
                var e = [];
                window.addEventListener("message", function(m) {
                    var h = m.source;
                    h !== window && null !== h || "process-tick" !== m.data || (m.stopPropagation(), 0 < e.length && e.shift()())
                }, !0);
                return function(m) {
                    e.push(m);
                    window.postMessage("process-tick", "*")
                }
            }
            return function(e) {
                setTimeout(e, 0)
            }
        }();
        e.title = "browser";
        e.browser = !0;
        e.env = {};
        e.argv = [];
        e.on = p;
        e.addListener = p;
        e.once = p;
        e.off = p;
        e.removeListener = p;
        e.removeAllListeners = p;
        e.emit = p;
        e.binding = function(e) {
            throw Error("process.binding is not supported");
        };
        e.cwd = function() {
            return "/"
        };
        e.chdir = function(e) {
            throw Error("process.chdir is not supported");
        }
    }, {}],
    64: [function(e, u, s) {
        function p() {}

        function n() {
            this.was = [0]
        }

        function m(a, d, f) {
            this.hufts = new Int32Array(3 * b);
            this.window = new Uint8Array(f);
            this.end = f;
            this.checkfn = d;
            this.mode = t;
            this.reset(a, null);
            this.index = this.table =
                this.left = 0;
            this.blens = null;
            this.bb = new Int32Array(1);
            this.tb = new Int32Array(1);
            this.codes = new h;
            this.check = this.write = this.read = this.bitb = this.bitk = this.last = 0;
            this.inftree = new v
        }

        function h() {}

        function v() {}

        function r(a, b, d, f, e) {
            if (0 != e) {
                if (!a) throw "Undef src";
                if (!d) throw "Undef dest";
                if (0 == b && e == a.length) d.set(a, f);
                else if (F) a = a.subarray(b, b + e), d.set(a, f);
                else if (1 == a.BYTES_PER_ELEMENT && 100 < e) a = new Uint8Array(a.buffer, a.byteOffset + b, e), d.set(a, f);
                else
                    for (var g = 0; g < e; ++g) d[f + g] = a[b + g]
            }
        }

        function z(c,
            b, d, e) {
            c = b ? d ? new Uint8Array(c, b, d) : new Uint8Array(c, b, c.byteLength - b) : new Uint8Array(c);
            d = new p;
            d.inflateInit(a, !0);
            d.next_in = c;
            d.next_in_index = 0;
            d.avail_in = c.length;
            c = [];
            for (var h = 0;;) {
                var n = new Uint8Array(32E3);
                d.next_out = n;
                d.next_out_index = 0;
                d.avail_out = n.length;
                var m = d.inflate(q);
                if (m != g && m != l && m != f) throw d.msg;
                if (0 != d.avail_out) {
                    var s = new Uint8Array(n.length - d.avail_out);
                    r(n, 0, s, 0, n.length - d.avail_out);
                    n = s
                }
                c.push(n);
                h += n.length;
                if (m == l || m == f) break
            }
            e && (e[0] = (b || 0) + d.next_in_index);
            if (1 == c.length) return c[0].buffer;
            b = new Uint8Array(h);
            for (d = e = 0; d < c.length; ++d) h = c[d], r(h, 0, b, e, h.length), e += h.length;
            return b.buffer
        }
        var a = 15,
            b = 1440,
            q = 0,
            g = 0,
            l = 1,
            f = -5,
            d = [0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535],
            t = 0,
            D = [96, 7, 256, 0, 8, 80, 0, 8, 16, 84, 8, 115, 82, 7, 31, 0, 8, 112, 0, 8, 48, 0, 9, 192, 80, 7, 10, 0, 8, 96, 0, 8, 32, 0, 9, 160, 0, 8, 0, 0, 8, 128, 0, 8, 64, 0, 9, 224, 80, 7, 6, 0, 8, 88, 0, 8, 24, 0, 9, 144, 83, 7, 59, 0, 8, 120, 0, 8, 56, 0, 9, 208, 81, 7, 17, 0, 8, 104, 0, 8, 40, 0, 9, 176, 0, 8, 8, 0, 8, 136, 0, 8, 72, 0, 9, 240, 80, 7, 4, 0, 8, 84, 0, 8, 20, 85, 8, 227, 83, 7, 43, 0, 8, 116, 0,
                8, 52, 0, 9, 200, 81, 7, 13, 0, 8, 100, 0, 8, 36, 0, 9, 168, 0, 8, 4, 0, 8, 132, 0, 8, 68, 0, 9, 232, 80, 7, 8, 0, 8, 92, 0, 8, 28, 0, 9, 152, 84, 7, 83, 0, 8, 124, 0, 8, 60, 0, 9, 216, 82, 7, 23, 0, 8, 108, 0, 8, 44, 0, 9, 184, 0, 8, 12, 0, 8, 140, 0, 8, 76, 0, 9, 248, 80, 7, 3, 0, 8, 82, 0, 8, 18, 85, 8, 163, 83, 7, 35, 0, 8, 114, 0, 8, 50, 0, 9, 196, 81, 7, 11, 0, 8, 98, 0, 8, 34, 0, 9, 164, 0, 8, 2, 0, 8, 130, 0, 8, 66, 0, 9, 228, 80, 7, 7, 0, 8, 90, 0, 8, 26, 0, 9, 148, 84, 7, 67, 0, 8, 122, 0, 8, 58, 0, 9, 212, 82, 7, 19, 0, 8, 106, 0, 8, 42, 0, 9, 180, 0, 8, 10, 0, 8, 138, 0, 8, 74, 0, 9, 244, 80, 7, 5, 0, 8, 86, 0, 8, 22, 192, 8, 0, 83, 7, 51, 0, 8, 118, 0, 8, 54, 0, 9, 204, 81, 7, 15,
                0, 8, 102, 0, 8, 38, 0, 9, 172, 0, 8, 6, 0, 8, 134, 0, 8, 70, 0, 9, 236, 80, 7, 9, 0, 8, 94, 0, 8, 30, 0, 9, 156, 84, 7, 99, 0, 8, 126, 0, 8, 62, 0, 9, 220, 82, 7, 27, 0, 8, 110, 0, 8, 46, 0, 9, 188, 0, 8, 14, 0, 8, 142, 0, 8, 78, 0, 9, 252, 96, 7, 256, 0, 8, 81, 0, 8, 17, 85, 8, 131, 82, 7, 31, 0, 8, 113, 0, 8, 49, 0, 9, 194, 80, 7, 10, 0, 8, 97, 0, 8, 33, 0, 9, 162, 0, 8, 1, 0, 8, 129, 0, 8, 65, 0, 9, 226, 80, 7, 6, 0, 8, 89, 0, 8, 25, 0, 9, 146, 83, 7, 59, 0, 8, 121, 0, 8, 57, 0, 9, 210, 81, 7, 17, 0, 8, 105, 0, 8, 41, 0, 9, 178, 0, 8, 9, 0, 8, 137, 0, 8, 73, 0, 9, 242, 80, 7, 4, 0, 8, 85, 0, 8, 21, 80, 8, 258, 83, 7, 43, 0, 8, 117, 0, 8, 53, 0, 9, 202, 81, 7, 13, 0, 8, 101, 0, 8, 37, 0,
                9, 170, 0, 8, 5, 0, 8, 133, 0, 8, 69, 0, 9, 234, 80, 7, 8, 0, 8, 93, 0, 8, 29, 0, 9, 154, 84, 7, 83, 0, 8, 125, 0, 8, 61, 0, 9, 218, 82, 7, 23, 0, 8, 109, 0, 8, 45, 0, 9, 186, 0, 8, 13, 0, 8, 141, 0, 8, 77, 0, 9, 250, 80, 7, 3, 0, 8, 83, 0, 8, 19, 85, 8, 195, 83, 7, 35, 0, 8, 115, 0, 8, 51, 0, 9, 198, 81, 7, 11, 0, 8, 99, 0, 8, 35, 0, 9, 166, 0, 8, 3, 0, 8, 131, 0, 8, 67, 0, 9, 230, 80, 7, 7, 0, 8, 91, 0, 8, 27, 0, 9, 150, 84, 7, 67, 0, 8, 123, 0, 8, 59, 0, 9, 214, 82, 7, 19, 0, 8, 107, 0, 8, 43, 0, 9, 182, 0, 8, 11, 0, 8, 139, 0, 8, 75, 0, 9, 246, 80, 7, 5, 0, 8, 87, 0, 8, 23, 192, 8, 0, 83, 7, 51, 0, 8, 119, 0, 8, 55, 0, 9, 206, 81, 7, 15, 0, 8, 103, 0, 8, 39, 0, 9, 174, 0, 8, 7, 0, 8, 135,
                0, 8, 71, 0, 9, 238, 80, 7, 9, 0, 8, 95, 0, 8, 31, 0, 9, 158, 84, 7, 99, 0, 8, 127, 0, 8, 63, 0, 9, 222, 82, 7, 27, 0, 8, 111, 0, 8, 47, 0, 9, 190, 0, 8, 15, 0, 8, 143, 0, 8, 79, 0, 9, 254, 96, 7, 256, 0, 8, 80, 0, 8, 16, 84, 8, 115, 82, 7, 31, 0, 8, 112, 0, 8, 48, 0, 9, 193, 80, 7, 10, 0, 8, 96, 0, 8, 32, 0, 9, 161, 0, 8, 0, 0, 8, 128, 0, 8, 64, 0, 9, 225, 80, 7, 6, 0, 8, 88, 0, 8, 24, 0, 9, 145, 83, 7, 59, 0, 8, 120, 0, 8, 56, 0, 9, 209, 81, 7, 17, 0, 8, 104, 0, 8, 40, 0, 9, 177, 0, 8, 8, 0, 8, 136, 0, 8, 72, 0, 9, 241, 80, 7, 4, 0, 8, 84, 0, 8, 20, 85, 8, 227, 83, 7, 43, 0, 8, 116, 0, 8, 52, 0, 9, 201, 81, 7, 13, 0, 8, 100, 0, 8, 36, 0, 9, 169, 0, 8, 4, 0, 8, 132, 0, 8, 68, 0, 9, 233, 80,
                7, 8, 0, 8, 92, 0, 8, 28, 0, 9, 153, 84, 7, 83, 0, 8, 124, 0, 8, 60, 0, 9, 217, 82, 7, 23, 0, 8, 108, 0, 8, 44, 0, 9, 185, 0, 8, 12, 0, 8, 140, 0, 8, 76, 0, 9, 249, 80, 7, 3, 0, 8, 82, 0, 8, 18, 85, 8, 163, 83, 7, 35, 0, 8, 114, 0, 8, 50, 0, 9, 197, 81, 7, 11, 0, 8, 98, 0, 8, 34, 0, 9, 165, 0, 8, 2, 0, 8, 130, 0, 8, 66, 0, 9, 229, 80, 7, 7, 0, 8, 90, 0, 8, 26, 0, 9, 149, 84, 7, 67, 0, 8, 122, 0, 8, 58, 0, 9, 213, 82, 7, 19, 0, 8, 106, 0, 8, 42, 0, 9, 181, 0, 8, 10, 0, 8, 138, 0, 8, 74, 0, 9, 245, 80, 7, 5, 0, 8, 86, 0, 8, 22, 192, 8, 0, 83, 7, 51, 0, 8, 118, 0, 8, 54, 0, 9, 205, 81, 7, 15, 0, 8, 102, 0, 8, 38, 0, 9, 173, 0, 8, 6, 0, 8, 134, 0, 8, 70, 0, 9, 237, 80, 7, 9, 0, 8, 94, 0, 8, 30, 0,
                9, 157, 84, 7, 99, 0, 8, 126, 0, 8, 62, 0, 9, 221, 82, 7, 27, 0, 8, 110, 0, 8, 46, 0, 9, 189, 0, 8, 14, 0, 8, 142, 0, 8, 78, 0, 9, 253, 96, 7, 256, 0, 8, 81, 0, 8, 17, 85, 8, 131, 82, 7, 31, 0, 8, 113, 0, 8, 49, 0, 9, 195, 80, 7, 10, 0, 8, 97, 0, 8, 33, 0, 9, 163, 0, 8, 1, 0, 8, 129, 0, 8, 65, 0, 9, 227, 80, 7, 6, 0, 8, 89, 0, 8, 25, 0, 9, 147, 83, 7, 59, 0, 8, 121, 0, 8, 57, 0, 9, 211, 81, 7, 17, 0, 8, 105, 0, 8, 41, 0, 9, 179, 0, 8, 9, 0, 8, 137, 0, 8, 73, 0, 9, 243, 80, 7, 4, 0, 8, 85, 0, 8, 21, 80, 8, 258, 83, 7, 43, 0, 8, 117, 0, 8, 53, 0, 9, 203, 81, 7, 13, 0, 8, 101, 0, 8, 37, 0, 9, 171, 0, 8, 5, 0, 8, 133, 0, 8, 69, 0, 9, 235, 80, 7, 8, 0, 8, 93, 0, 8, 29, 0, 9, 155, 84, 7, 83, 0, 8,
                125, 0, 8, 61, 0, 9, 219, 82, 7, 23, 0, 8, 109, 0, 8, 45, 0, 9, 187, 0, 8, 13, 0, 8, 141, 0, 8, 77, 0, 9, 251, 80, 7, 3, 0, 8, 83, 0, 8, 19, 85, 8, 195, 83, 7, 35, 0, 8, 115, 0, 8, 51, 0, 9, 199, 81, 7, 11, 0, 8, 99, 0, 8, 35, 0, 9, 167, 0, 8, 3, 0, 8, 131, 0, 8, 67, 0, 9, 231, 80, 7, 7, 0, 8, 91, 0, 8, 27, 0, 9, 151, 84, 7, 67, 0, 8, 123, 0, 8, 59, 0, 9, 215, 82, 7, 19, 0, 8, 107, 0, 8, 43, 0, 9, 183, 0, 8, 11, 0, 8, 139, 0, 8, 75, 0, 9, 247, 80, 7, 5, 0, 8, 87, 0, 8, 23, 192, 8, 0, 83, 7, 51, 0, 8, 119, 0, 8, 55, 0, 9, 207, 81, 7, 15, 0, 8, 103, 0, 8, 39, 0, 9, 175, 0, 8, 7, 0, 8, 135, 0, 8, 71, 0, 9, 239, 80, 7, 9, 0, 8, 95, 0, 8, 31, 0, 9, 159, 84, 7, 99, 0, 8, 127, 0, 8, 63, 0, 9, 223,
                82, 7, 27, 0, 8, 111, 0, 8, 47, 0, 9, 191, 0, 8, 15, 0, 8, 143, 0, 8, 79, 0, 9, 255
            ],
            P = [80, 5, 1, 87, 5, 257, 83, 5, 17, 91, 5, 4097, 81, 5, 5, 89, 5, 1025, 85, 5, 65, 93, 5, 16385, 80, 5, 3, 88, 5, 513, 84, 5, 33, 92, 5, 8193, 82, 5, 9, 90, 5, 2049, 86, 5, 129, 192, 5, 24577, 80, 5, 2, 87, 5, 385, 83, 5, 25, 91, 5, 6145, 81, 5, 7, 89, 5, 1537, 85, 5, 97, 93, 5, 24577, 80, 5, 4, 88, 5, 769, 84, 5, 49, 92, 5, 12289, 82, 5, 13, 90, 5, 3073, 86, 5, 193, 192, 5, 24577],
            G = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0],
            L = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,
                5, 5, 0, 112, 112
            ],
            N = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577],
            Q = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
        p.prototype.inflateInit = function(c, b) {
            c || (c = a);
            b && (b = !1);
            this.istate = new n;
            return this.istate.inflateInit(this, b ? -c : c)
        };
        p.prototype.inflate = function(a) {
            return null == this.istate ? -2 : this.istate.inflate(this, a)
        };
        p.prototype.inflateEnd = function() {
            if (null == this.istate) return -2;
            var a = istate.inflateEnd(this);
            this.istate = null;
            return a
        };
        p.prototype.inflateSync = function() {
            return istate.inflateSync(this)
        };
        p.prototype.inflateSetDictionary = function(a, b) {
            return istate.inflateSetDictionary(this, a, b)
        };
        n.prototype.inflateReset = function(a) {
            if (null == a || null == a.istate) return -2;
            a.total_in = a.total_out = 0;
            a.msg = null;
            a.istate.mode = 0 != a.istate.nowrap ? 7 : 0;
            a.istate.blocks.reset(a, null);
            return g
        };
        n.prototype.inflateEnd = function(a) {
            null != this.blocks && this.blocks.free(a);
            this.blocks = null;
            return g
        };
        n.prototype.inflateInit = function(a,
            b) {
            this.blocks = a.msg = null;
            nowrap = 0;
            0 > b && (b = -b, nowrap = 1);
            if (8 > b || 15 < b) return this.inflateEnd(a), -2;
            this.wbits = b;
            a.istate.blocks = new m(a, 0 != a.istate.nowrap ? null : this, 1 << b);
            this.inflateReset(a);
            return g
        };
        n.prototype.inflate = function(a, b) {
            var d, e;
            if (null == a || null == a.istate || null == a.next_in) return -2;
            b = 4 == b ? f : g;
            for (d = f;;) switch (a.istate.mode) {
                case 0:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    if (8 != ((a.istate.method = a.next_in[a.next_in_index++]) & 15)) {
                        a.istate.mode = 13;
                        a.msg = "unknown compression method";
                        a.istate.marker = 5;
                        break
                    }
                    if ((a.istate.method >> 4) + 8 > a.istate.wbits) {
                        a.istate.mode = 13;
                        a.msg = "invalid window size";
                        a.istate.marker = 5;
                        break
                    }
                    a.istate.mode = 1;
                case 1:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    e = a.next_in[a.next_in_index++] & 255;
                    if (0 != ((a.istate.method << 8) + e) % 31) {
                        a.istate.mode = 13;
                        a.msg = "incorrect header check";
                        a.istate.marker = 5;
                        break
                    }
                    if (0 == (e & 32)) {
                        a.istate.mode = 7;
                        break
                    }
                    a.istate.mode = 2;
                case 2:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need = (a.next_in[a.next_in_index++] &
                        255) << 24 & 4278190080;
                    a.istate.mode = 3;
                case 3:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need += (a.next_in[a.next_in_index++] & 255) << 16 & 16711680;
                    a.istate.mode = 4;
                case 4:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need += (a.next_in[a.next_in_index++] & 255) << 8 & 65280;
                    a.istate.mode = 5;
                case 5:
                    if (0 == a.avail_in) return d;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need += a.next_in[a.next_in_index++] & 255;
                    a.adler = a.istate.need;
                    a.istate.mode = 6;
                    return 2;
                case 6:
                    return a.istate.mode = 13,
                        a.msg = "need dictionary", a.istate.marker = 0, -2;
                case 7:
                    d = a.istate.blocks.proc(a, d);
                    if (-3 == d) {
                        a.istate.mode = 13;
                        a.istate.marker = 0;
                        break
                    }
                    d == g && (d = b);
                    if (d != l) return d;
                    d = b;
                    a.istate.blocks.reset(a, a.istate.was);
                    if (0 != a.istate.nowrap) {
                        a.istate.mode = 12;
                        break
                    }
                    a.istate.mode = 8;
                case 8:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need = (a.next_in[a.next_in_index++] & 255) << 24 & 4278190080;
                    a.istate.mode = 9;
                case 9:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need += (a.next_in[a.next_in_index++] &
                        255) << 16 & 16711680;
                    a.istate.mode = 10;
                case 10:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need += (a.next_in[a.next_in_index++] & 255) << 8 & 65280;
                    a.istate.mode = 11;
                case 11:
                    if (0 == a.avail_in) return d;
                    d = b;
                    a.avail_in--;
                    a.total_in++;
                    a.istate.need += a.next_in[a.next_in_index++] & 255;
                    if (a.istate.was[0] != a.istate.need) {
                        a.istate.mode = 13;
                        a.msg = "incorrect data check";
                        a.istate.marker = 5;
                        break
                    }
                    a.istate.mode = 12;
                case 12:
                    return l;
                case 13:
                    return -3;
                default:
                    return -2
            }
        };
        n.prototype.inflateSetDictionary = function(a,
            b, d) {
            var e = 0,
                f = d;
            if (null == a || null == a.istate || 6 != a.istate.mode) return -2;
            if (a._adler.adler32(1, b, 0, d) != a.adler) return -3;
            a.adler = a._adler.adler32(0, null, 0, 0);
            f >= 1 << a.istate.wbits && (f = (1 << a.istate.wbits) - 1, e = d - f);
            a.istate.blocks.set_dictionary(b, e, f);
            a.istate.mode = 7;
            return g
        };
        var y = [0, 0, 255, 255];
        n.prototype.inflateSync = function(a) {
            var b, d, e;
            if (null == a || null == a.istate) return -2;
            13 != a.istate.mode && (a.istate.mode = 13, a.istate.marker = 0);
            if (0 == (b = a.avail_in)) return f;
            d = a.next_in_index;
            for (e = a.istate.marker; 0 !=
                b && 4 > e;) a.next_in[d] == y[e] ? e++ : e = 0 != a.next_in[d] ? 0 : 4 - e, d++, b--;
            a.total_in += d - a.next_in_index;
            a.next_in_index = d;
            a.avail_in = b;
            a.istate.marker = e;
            if (4 != e) return -3;
            b = a.total_in;
            d = a.total_out;
            this.inflateReset(a);
            a.total_in = b;
            a.total_out = d;
            a.istate.mode = 7;
            return g
        };
        n.prototype.inflateSyncPoint = function(a) {
            return null == a || null == a.istate || null == a.istate.blocks ? -2 : a.istate.blocks.sync_point()
        };
        var H = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
        m.prototype.reset = function(a, b) {
            b && (b[0] = this.check);
            6 == this.mode &&
                this.codes.free(a);
            this.mode = t;
            this.read = this.write = this.bitb = this.bitk = 0;
            this.checkfn && (a.adler = this.check = a._adler.adler32(0, null, 0, 0))
        };
        m.prototype.proc = function(a, b) {
            var e, f, h, n, m, p, k;
            n = a.next_in_index;
            m = a.avail_in;
            f = this.bitb;
            h = this.bitk;
            p = this.write;
            for (k = p < this.read ? this.read - p - 1 : this.end - p;;) switch (this.mode) {
                case t:
                    for (; 3 > h;) {
                        if (0 != m) b = g;
                        else return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                        m--;
                        f |= (a.next_in[n++] &
                            255) << h;
                        h += 8
                    }
                    e = f & 7;
                    this.last = e & 1;
                    switch (e >>> 1) {
                        case 0:
                            f >>>= 3;
                            h -= 3;
                            e = h & 7;
                            f >>>= e;
                            h -= e;
                            this.mode = 1;
                            break;
                        case 1:
                            var q = new Int32Array(1),
                                s = new Int32Array(1),
                                u = [],
                                v = [];
                            e = s;
                            var y = u,
                                z = v;
                            q[0] = 9;
                            e[0] = 5;
                            y[0] = D;
                            z[0] = P;
                            this.codes.init(q[0], s[0], u[0], 0, v[0], 0, a);
                            f >>>= 3;
                            h -= 3;
                            this.mode = 6;
                            break;
                        case 2:
                            f >>>= 3;
                            h -= 3;
                            this.mode = 3;
                            break;
                        case 3:
                            return f >>>= 3, h -= 3, this.mode = 13, a.msg = "invalid block type", b = -3, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a,
                                b)
                    }
                    break;
                case 1:
                    for (; 32 > h;) {
                        if (0 != m) b = g;
                        else return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                        m--;
                        f |= (a.next_in[n++] & 255) << h;
                        h += 8
                    }
                    if ((~f >>> 16 & 65535) != (f & 65535)) return this.mode = 13, a.msg = "invalid stored block lengths", b = -3, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                    this.left = f & 65535;
                    f = h = 0;
                    this.mode = 0 != this.left ? 2 : 0 != this.last ? 7 : t;
                    break;
                case 2:
                    if (0 ==
                        m) return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, write = p, this.inflate_flush(a, b);
                    if (0 == k && (p == end && 0 != read && (p = 0, k = p < this.read ? this.read - p - 1 : this.end - p), 0 == k && (this.write = p, b = this.inflate_flush(a, b), p = this.write, k = p < this.read ? this.read - p - 1 : this.end - p, p == this.end && 0 != this.read && (p = 0, k = p < this.read ? this.read - p - 1 : this.end - p), 0 == k))) return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a,
                        b);
                    b = g;
                    e = this.left;
                    e > m && (e = m);
                    e > k && (e = k);
                    r(a.next_in, n, this.window, p, e);
                    n += e;
                    m -= e;
                    p += e;
                    k -= e;
                    if (0 != (this.left -= e)) break;
                    this.mode = 0 != this.last ? 7 : t;
                    break;
                case 3:
                    for (; 14 > h;) {
                        if (0 != m) b = g;
                        else return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                        m--;
                        f |= (a.next_in[n++] & 255) << h;
                        h += 8
                    }
                    this.table = e = f & 16383;
                    if (29 < (e & 31) || 29 < (e >> 5 & 31)) return this.mode = 9, a.msg = "too many length or distance symbols", b = -3, this.bitb = f, this.bitk = h, a.avail_in =
                        m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                    e = 258 + (e & 31) + (e >> 5 & 31);
                    if (null == this.blens || this.blens.length < e) this.blens = new Int32Array(e);
                    else
                        for (k = 0; k < e; k++) this.blens[k] = 0;
                    f >>>= 14;
                    h -= 14;
                    this.index = 0;
                    mode = 4;
                case 4:
                    for (; this.index < 4 + (this.table >>> 10);) {
                        for (; 3 > h;) {
                            if (0 != m) b = g;
                            else return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                            m--;
                            f |= (a.next_in[n++] & 255) << h;
                            h += 8
                        }
                        this.blens[H[this.index++]] =
                            f & 7;
                        f >>>= 3;
                        h -= 3
                    }
                    for (; 19 > this.index;) this.blens[H[this.index++]] = 0;
                    this.bb[0] = 7;
                    e = this.inftree.inflate_trees_bits(this.blens, this.bb, this.tb, this.hufts, a);
                    if (e != g) return b = e, -3 == b && (this.blens = null, this.mode = 9), this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, write = p, this.inflate_flush(a, b);
                    this.index = 0;
                    this.mode = 5;
                case 5:
                    for (;;) {
                        e = this.table;
                        if (!(this.index < 258 + (e & 31) + (e >> 5 & 31))) break;
                        for (e = this.bb[0]; h < e;) {
                            if (0 != m) b = g;
                            else return this.bitb = f, this.bitk = h, a.avail_in =
                                m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                            m--;
                            f |= (a.next_in[n++] & 255) << h;
                            h += 8
                        }
                        e = this.hufts[3 * (this.tb[0] + (f & d[e])) + 1];
                        s = this.hufts[3 * (this.tb[0] + (f & d[e])) + 2];
                        if (16 > s) f >>>= e, h -= e, this.blens[this.index++] = s;
                        else {
                            k = 18 == s ? 7 : s - 14;
                            for (q = 18 == s ? 11 : 3; h < e + k;) {
                                if (0 != m) b = g;
                                else return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                                m--;
                                f |= (a.next_in[n++] & 255) << h;
                                h += 8
                            }
                            f >>>= e;
                            h -= e;
                            q += f & d[k];
                            f >>>=
                                k;
                            h -= k;
                            k = this.index;
                            e = this.table;
                            if (k + q > 258 + (e & 31) + (e >> 5 & 31) || 16 == s && 1 > k) return this.blens = null, this.mode = 9, a.msg = "invalid bit length repeat", b = -3, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                            s = 16 == s ? this.blens[k - 1] : 0;
                            do this.blens[k++] = s; while (0 != --q);
                            this.index = k
                        }
                    }
                    this.tb[0] = -1;
                    q = new Int32Array(1);
                    s = new Int32Array(1);
                    u = new Int32Array(1);
                    v = new Int32Array(1);
                    q[0] = 9;
                    s[0] = 6;
                    e = this.table;
                    e = this.inftree.inflate_trees_dynamic(257 +
                        (e & 31), 1 + (e >> 5 & 31), this.blens, q, s, u, v, this.hufts, a);
                    if (e != g) return -3 == e && (this.blens = null, this.mode = 13), b = e, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                    this.codes.init(q[0], s[0], this.hufts, u[0], this.hufts, v[0], a);
                    this.mode = 6;
                case 6:
                    this.bitb = f;
                    this.bitk = h;
                    a.avail_in = m;
                    a.total_in += n - a.next_in_index;
                    a.next_in_index = n;
                    this.write = p;
                    if ((b = this.codes.proc(this, a, b)) != l) return this.inflate_flush(a, b);
                    b = g;
                    this.codes.free(a);
                    n = a.next_in_index;
                    m = a.avail_in;
                    f = this.bitb;
                    h = this.bitk;
                    p = this.write;
                    k = p < this.read ? this.read - p - 1 : this.end - p;
                    if (0 == this.last) {
                        this.mode = t;
                        break
                    }
                    this.mode = 7;
                case 7:
                    this.write = p;
                    b = this.inflate_flush(a, b);
                    p = this.write;
                    if (this.read != this.write) return this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                    mode = 12;
                case 8:
                    return b = l, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                case 9:
                    return b = -3, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b);
                default:
                    return b = -2, this.bitb = f, this.bitk = h, a.avail_in = m, a.total_in += n - a.next_in_index, a.next_in_index = n, this.write = p, this.inflate_flush(a, b)
            }
        };
        m.prototype.free = function(a) {
            this.reset(a, null);
            this.hufts = this.window = null
        };
        m.prototype.set_dictionary = function(a, b, d) {
            r(a, b, window, 0, d);
            this.read = this.write = d
        };
        m.prototype.sync_point = function() {
            return 1 == this.mode
        };
        m.prototype.inflate_flush =
            function(a, b) {
                var d, e, h;
                e = a.next_out_index;
                h = this.read;
                d = (h <= this.write ? this.write : this.end) - h;
                d > a.avail_out && (d = a.avail_out);
                0 != d && b == f && (b = g);
                a.avail_out -= d;
                a.total_out += d;
                null != this.checkfn && (a.adler = this.check = a._adler.adler32(this.check, this.window, h, d));
                r(this.window, h, a.next_out, e, d);
                e += d;
                h += d;
                h == this.end && (h = 0, this.write == this.end && (this.write = 0), d = this.write - h, d > a.avail_out && (d = a.avail_out), 0 != d && b == f && (b = g), a.avail_out -= d, a.total_out += d, null != this.checkfn && (a.adler = this.check = a._adler.adler32(this.check,
                    this.window, h, d)), r(this.window, h, a.next_out, e, d), e += d, h += d);
                a.next_out_index = e;
                this.read = h;
                return b
            };
        h.prototype.init = function(a, b, d, e, f, g, h) {
            this.mode = 0;
            this.lbits = a;
            this.dbits = b;
            this.ltree = d;
            this.ltree_index = e;
            this.dtree = f;
            this.dtree_index = g;
            this.tree = null
        };
        h.prototype.proc = function(a, b, e) {
            var f, h, n = 0,
                m = 0,
                p = 0,
                k, q, r, p = b.next_in_index;
            k = b.avail_in;
            n = a.bitb;
            m = a.bitk;
            q = a.write;
            for (r = q < a.read ? a.read - q - 1 : a.end - q;;) switch (this.mode) {
                case 0:
                    if (258 <= r && 10 <= k && (a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in +=
                            p - b.next_in_index, b.next_in_index = p, a.write = q, e = this.inflate_fast(this.lbits, this.dbits, this.ltree, this.ltree_index, this.dtree, this.dtree_index, a, b), p = b.next_in_index, k = b.avail_in, n = a.bitb, m = a.bitk, q = a.write, r = q < a.read ? a.read - q - 1 : a.end - q, e != g)) {
                        this.mode = e == l ? 7 : 9;
                        break
                    }
                    this.need = this.lbits;
                    this.tree = this.ltree;
                    this.tree_index = this.ltree_index;
                    this.mode = 1;
                case 1:
                    for (f = this.need; m < f;) {
                        if (0 != k) e = g;
                        else return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b,
                            e);
                        k--;
                        n |= (b.next_in[p++] & 255) << m;
                        m += 8
                    }
                    f = 3 * (this.tree_index + (n & d[f]));
                    n >>>= this.tree[f + 1];
                    m -= this.tree[f + 1];
                    h = this.tree[f];
                    if (0 == h) {
                        this.lit = this.tree[f + 2];
                        this.mode = 6;
                        break
                    }
                    if (0 != (h & 16)) {
                        this.get = h & 15;
                        this.len = this.tree[f + 2];
                        this.mode = 2;
                        break
                    }
                    if (0 == (h & 64)) {
                        this.need = h;
                        this.tree_index = f / 3 + this.tree[f + 2];
                        break
                    }
                    if (0 != (h & 32)) {
                        this.mode = 7;
                        break
                    }
                    this.mode = 9;
                    b.msg = "invalid literal/length code";
                    e = -3;
                    a.bitb = n;
                    a.bitk = m;
                    b.avail_in = k;
                    b.total_in += p - b.next_in_index;
                    b.next_in_index = p;
                    a.write = q;
                    return a.inflate_flush(b,
                        e);
                case 2:
                    for (f = this.get; m < f;) {
                        if (0 != k) e = g;
                        else return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                        k--;
                        n |= (b.next_in[p++] & 255) << m;
                        m += 8
                    }
                    this.len += n & d[f];
                    n >>= f;
                    m -= f;
                    this.need = this.dbits;
                    this.tree = this.dtree;
                    this.tree_index = this.dtree_index;
                    this.mode = 3;
                case 3:
                    for (f = this.need; m < f;) {
                        if (0 != k) e = g;
                        else return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                        k--;
                        n |= (b.next_in[p++] &
                            255) << m;
                        m += 8
                    }
                    f = 3 * (this.tree_index + (n & d[f]));
                    n >>= this.tree[f + 1];
                    m -= this.tree[f + 1];
                    h = this.tree[f];
                    if (0 != (h & 16)) {
                        this.get = h & 15;
                        this.dist = this.tree[f + 2];
                        this.mode = 4;
                        break
                    }
                    if (0 == (h & 64)) {
                        this.need = h;
                        this.tree_index = f / 3 + this.tree[f + 2];
                        break
                    }
                    this.mode = 9;
                    b.msg = "invalid distance code";
                    e = -3;
                    a.bitb = n;
                    a.bitk = m;
                    b.avail_in = k;
                    b.total_in += p - b.next_in_index;
                    b.next_in_index = p;
                    a.write = q;
                    return a.inflate_flush(b, e);
                case 4:
                    for (f = this.get; m < f;) {
                        if (0 != k) e = g;
                        else return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index,
                            b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                        k--;
                        n |= (b.next_in[p++] & 255) << m;
                        m += 8
                    }
                    this.dist += n & d[f];
                    n >>= f;
                    m -= f;
                    this.mode = 5;
                case 5:
                    for (f = q - this.dist; 0 > f;) f += a.end;
                    for (; 0 != this.len;) {
                        if (0 == r && (q == a.end && 0 != a.read && (q = 0, r = q < a.read ? a.read - q - 1 : a.end - q), 0 == r && (a.write = q, e = a.inflate_flush(b, e), q = a.write, r = q < a.read ? a.read - q - 1 : a.end - q, q == a.end && 0 != a.read && (q = 0, r = q < a.read ? a.read - q - 1 : a.end - q), 0 == r))) return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b,
                            e);
                        a.window[q++] = a.window[f++];
                        r--;
                        f == a.end && (f = 0);
                        this.len--
                    }
                    this.mode = 0;
                    break;
                case 6:
                    if (0 == r && (q == a.end && 0 != a.read && (q = 0, r = q < a.read ? a.read - q - 1 : a.end - q), 0 == r && (a.write = q, e = a.inflate_flush(b, e), q = a.write, r = q < a.read ? a.read - q - 1 : a.end - q, q == a.end && 0 != a.read && (q = 0, r = q < a.read ? a.read - q - 1 : a.end - q), 0 == r))) return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                    e = g;
                    a.window[q++] = this.lit;
                    r--;
                    this.mode = 0;
                    break;
                case 7:
                    7 < m && (m -= 8, k++, p--);
                    a.write = q;
                    e = a.inflate_flush(b, e);
                    q = a.write;
                    if (a.read != a.write) return a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                    this.mode = 8;
                case 8:
                    return e = l, a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                case 9:
                    return e = -3, a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in += p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e);
                default:
                    return e = -2, a.bitb = n, a.bitk = m, b.avail_in = k, b.total_in +=
                        p - b.next_in_index, b.next_in_index = p, a.write = q, a.inflate_flush(b, e)
            }
        };
        h.prototype.free = function(a) {};
        h.prototype.inflate_fast = function(a, b, e, f, h, n, m, p) {
            var k, q, s, u, t, v, y, z, F, H, D, G;
            v = p.next_in_index;
            y = p.avail_in;
            u = m.bitb;
            t = m.bitk;
            z = m.write;
            F = z < m.read ? m.read - z - 1 : m.end - z;
            a = d[a];
            H = d[b];
            do {
                for (; 20 > t;) y--, u |= (p.next_in[v++] & 255) << t, t += 8;
                k = u & a;
                q = e;
                s = f;
                G = 3 * (s + k);
                if (0 == (b = q[G])) u >>= q[G + 1], t -= q[G + 1], m.window[z++] = q[G + 2], F--;
                else {
                    do {
                        u >>= q[G + 1];
                        t -= q[G + 1];
                        if (0 != (b & 16)) {
                            b &= 15;
                            D = q[G + 2] + (u & d[b]);
                            u >>= b;
                            for (t -= b; 15 >
                                t;) y--, u |= (p.next_in[v++] & 255) << t, t += 8;
                            k = u & H;
                            q = h;
                            s = n;
                            G = 3 * (s + k);
                            b = q[G];
                            do
                                if (u >>= q[G + 1], t -= q[G + 1], 0 != (b & 16)) {
                                    for (b &= 15; t < b;) y--, u |= (p.next_in[v++] & 255) << t, t += 8;
                                    k = q[G + 2] + (u & d[b]);
                                    u >>= b;
                                    t -= b;
                                    F -= D;
                                    if (z >= k) k = z - k, m.window[z++] = m.window[k++], m.window[z++] = m.window[k++], D -= 2;
                                    else {
                                        k = z - k;
                                        do k += m.end; while (0 > k);
                                        b = m.end - k;
                                        if (D > b) {
                                            D -= b;
                                            if (0 < z - k && b > z - k) {
                                                do m.window[z++] = m.window[k++]; while (0 != --b)
                                            } else r(m.window, k, m.window, z, b), z += b;
                                            k = 0
                                        }
                                    }
                                    do m.window[z++] = m.window[k++]; while (0 != --D);
                                    break
                                } else if (0 == (b & 64)) k +=
                                q[G + 2], k += u & d[b], G = 3 * (s + k), b = q[G];
                            else return p.msg = "invalid distance code", D = p.avail_in - y, D = t >> 3 < D ? t >> 3 : D, y += D, v -= D, t -= D << 3, m.bitb = u, m.bitk = t, p.avail_in = y, p.total_in += v - p.next_in_index, p.next_in_index = v, m.write = z, -3;
                            while (1);
                            break
                        }
                        if (0 == (b & 64)) {
                            if (k += q[G + 2], k += u & d[b], G = 3 * (s + k), 0 == (b = q[G])) {
                                u >>= q[G + 1];
                                t -= q[G + 1];
                                m.window[z++] = q[G + 2];
                                F--;
                                break
                            }
                        } else {
                            if (0 != (b & 32)) return D = p.avail_in - y, D = t >> 3 < D ? t >> 3 : D, y += D, v -= D, t -= D << 3, m.bitb = u, m.bitk = t, p.avail_in = y, p.total_in += v - p.next_in_index, p.next_in_index = v, m.write =
                                z, l;
                            p.msg = "invalid literal/length code";
                            D = p.avail_in - y;
                            D = t >> 3 < D ? t >> 3 : D;
                            y += D;
                            v -= D;
                            t -= D << 3;
                            m.bitb = u;
                            m.bitk = t;
                            p.avail_in = y;
                            p.total_in += v - p.next_in_index;
                            p.next_in_index = v;
                            m.write = z;
                            return -3
                        }
                    } while (1)
                }
            } while (258 <= F && 10 <= y);
            D = p.avail_in - y;
            D = t >> 3 < D ? t >> 3 : D;
            v -= D;
            m.bitb = u;
            m.bitk = t - (D << 3);
            p.avail_in = y + D;
            p.total_in += v - p.next_in_index;
            p.next_in_index = v;
            m.write = z;
            return g
        };
        v.prototype.huft_build = function(a, d, e, h, l, m, n, p, k, q, s) {
            var u, t, v, y, z, D, F, G, H;
            D = 0;
            t = e;
            do this.c[a[d + D]]++, D++, t--; while (0 != t);
            if (this.c[0] ==
                e) return n[0] = -1, p[0] = 0, g;
            z = p[0];
            for (v = 1; 15 >= v && 0 == this.c[v]; v++);
            y = v;
            z < v && (z = v);
            for (t = 15; 0 != t && 0 == this.c[t]; t--);
            q = t;
            z > t && (z = t);
            p[0] = z;
            for (p = 1 << v; v < t; v++, p <<= 1)
                if (0 > (p -= this.c[v])) return -3;
            if (0 > (p -= this.c[t])) return -3;
            this.c[t] += p;
            this.x[1] = v = 0;
            D = 1;
            for (F = 2; 0 != --t;) this.x[F] = v += this.c[D], F++, D++;
            D = t = 0;
            do 0 != (v = a[d + D]) && (this.v[this.x[v]++] = t), D++; while (++t < e);
            e = this.x[q];
            D = this.x[0] = t = 0;
            d = -1;
            G = -z;
            for (H = F = this.u[0] = 0; y <= q; y++)
                for (a = this.c[y]; 0 != a--;) {
                    for (; y > G + z;) {
                        d++;
                        G += z;
                        H = q - G;
                        H = H > z ? z : H;
                        if ((u = 1 <<
                                (v = y - G)) > a + 1 && (u -= a + 1, F = y, v < H))
                            for (; ++v < H && !((u <<= 1) <= this.c[++F]);) u -= this.c[F];
                        H = 1 << v;
                        if (this.hn[0] + H > b) return -3;
                        this.u[d] = F = this.hn[0];
                        this.hn[0] += H;
                        0 != d ? (this.x[d] = t, this.r[0] = v, this.r[1] = z, v = t >>> G - z, this.r[2] = F - this.u[d - 1] - v, r(this.r, 0, k, 3 * (this.u[d - 1] + v), 3)) : n[0] = F
                    }
                    this.r[1] = y - G;
                    D >= e ? this.r[0] = 192 : s[D] < h ? (this.r[0] = 256 > this.v[D] ? 0 : 96, this.r[2] = this.v[D++]) : (this.r[0] = m[this.v[D] - h] + 16 + 64, this.r[2] = l[this.v[D++] - h]);
                    u = 1 << y - G;
                    for (v = t >>> G; v < H; v += u) r(this.r, 0, k, 3 * (F + v), 3);
                    for (v = 1 << y - 1; 0 != (t & v); v >>>=
                        1) t ^= v;
                    t ^= v;
                    for (v = (1 << G) - 1;
                        (t & v) != this.x[d];) d--, G -= z, v = (1 << G) - 1
                }
            return 0 != p && 1 != q ? f : g
        };
        v.prototype.inflate_trees_bits = function(a, b, d, e, g) {
            this.initWorkArea(19);
            this.hn[0] = 0;
            a = this.huft_build(a, 0, 19, 19, null, null, d, b, e, this.hn, this.v);
            if (-3 == a) g.msg = "oversubscribed dynamic bit lengths tree";
            else if (a == f || 0 == b[0]) g.msg = "incomplete dynamic bit lengths tree", a = -3;
            return a
        };
        v.prototype.inflate_trees_dynamic = function(a, b, d, e, h, l, m, n, k) {
            this.initWorkArea(288);
            this.hn[0] = 0;
            l = this.huft_build(d, 0, a, 257, G, L,
                l, e, n, this.hn, this.v);
            if (l != g || 0 == e[0]) return -3 == l ? k.msg = "oversubscribed literal/length tree" : -4 != l && (k.msg = "incomplete literal/length tree", l = -3), l;
            this.initWorkArea(288);
            l = this.huft_build(d, a, b, 0, N, Q, m, h, n, this.hn, this.v);
            return l != g || 0 == h[0] && 257 < a ? (-3 == l ? k.msg = "oversubscribed distance tree" : l == f ? (k.msg = "incomplete distance tree", l = -3) : -4 != l && (k.msg = "empty distance tree with lengths", l = -3), l) : g
        };
        v.prototype.initWorkArea = function(a) {
            null == this.hn && (this.hn = new Int32Array(1), this.v = new Int32Array(a),
                this.c = new Int32Array(16), this.r = new Int32Array(3), this.u = new Int32Array(15), this.x = new Int32Array(16));
            this.v.length < a && (this.v = new Int32Array(a));
            for (var b = 0; b < a; b++) this.v[b] = 0;
            for (b = 0; 16 > b; b++) this.c[b] = 0;
            for (b = 0; 3 > b; b++) this.r[b] = 0;
            r(this.c, 0, this.u, 0, 15);
            r(this.c, 0, this.x, 0, 16)
        };
        var F = "function" === typeof(new Uint8Array(1)).subarray;
        "undefined" !== typeof u && (u.exports = {
            inflateBuffer: z,
            arrayCopy: r
        })
    }, {}]
}, {}, [17]);
