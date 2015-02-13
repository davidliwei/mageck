/*global $, SF*/
/*jshint evil:true */

(function() {
    /******************************************************************
     Third party functions - included so webtracker can function
     standalone, without any other JS loaded.
     *****************************************************************/
    /*
        Find the value of a particular cookie
    */
    function munch(n) {
        var c = document.cookie;
        if (c) {
            var i = c.indexOf(n + '=');
            if (i > -1){
                var j = c.indexOf(';', i);
                return c.substring(i + n.length + 1, j < 0 ? c.length : j);
            }
        }
    }
    /*
        Developed by Robert Nyman, http://www.robertnyman.com
        Code/licensing: http://code.google.com/p/getelementsbyclassname/
    */
    var getElementsByClassName = function(b,a,c){if(document.getElementsByClassName){getElementsByClassName=function(j,m,h){h=h||document;var d=h.getElementsByClassName(j),l=(m)?new RegExp("\\b"+m+"\\b","i"):null,e=[],g;for(var f=0,k=d.length;f<k;f+=1){g=d[f];if(!l||l.test(g.nodeName)){e.push(g);}}return e;};}else{if(document.evaluate){getElementsByClassName=function(o,r,n){r=r||"*";n=n||document;var g=o.split(" "),p="",l="http://www.w3.org/1999/xhtml",q=(document.documentElement.namespaceURI===l)?l:null,h=[],d,f;for(var i=0,k=g.length;i<k;i+=1){p+="[contains(concat(' ', @class, ' '), ' "+g[i]+" ')]";}try{d=document.evaluate(".//"+r+p,n,q,0,null);}catch(m){d=document.evaluate(".//"+r+p,n,null,0,null);}while((f=d.iterateNext())){h.push(f);}return h;};}else{getElementsByClassName=function(r,u,q){u=u||"*";q=q||document;var h=r.split(" "),t=[],d=(u==="*"&&q.all)?q.all:q.getElementsByTagName(u),p,j=[],o;for(var i=0,e=h.length;i<e;i+=1){t.push(new RegExp("(^|\\s)"+h[i]+"(\\s|$)"));}for(var g=0,s=d.length;g<s;g+=1){p=d[g];o=false;for(var f=0,n=t.length;f<n;f+=1){o=t[f].test(p.className);if(!o){break;}}if(o){j.push(p);}}return j;};}}return getElementsByClassName(b,a,c);};
    /*
        From http://stackoverflow.com/questions/647259/javascript-query-string
    */
    function parseQueryString() {
        var result = {}, queryString = location.search.substring(1),
            re = /([^&=]+)=([^&]*)/g, m;
        while ((m = re.exec(queryString))) {
            result[decodeURIComponent(m[1])] = decodeURIComponent(m[2]);
        }
        return result;
    }
    /*
     * From http://stackoverflow.com/questions/4821627/building-a-google-analytics-domain-hash
     */
    function hash(d){var a=1,c=0,h,o;if(d){a=0;for(h=d.length-1;h>=0;h--){o=d.charCodeAt(h);a=(a<<6&268435455)+o+(o<<14);c=a&266338304;a=c!==0?a^c>>21:a;}} return a;} // jshint ignore:line

    /******************************************************************
     SF webtracking code
     *****************************************************************/
    var doc = document, body = doc.body,
        ads = getElementsByClassName('ad'),
        qs = parseQueryString(),
        url = '/log/webtracker/',
        testString = munch('__utmx') || munch('switchboard.test'),
        domainHash = munch('__utmc') || hash(window.location.hostname),
        data;
    // don't track error pages
    if (body.id === 'error-content') {
        return;
    }
    // parse json stored in meta element
    if (window.JSON && JSON.parse) {
        data = JSON.parse(doc.getElementById('webtracker').content);
    } else {
        data = eval('(' + doc.getElementById('webtracker').content + ')');
    }
    // grab url
    data.url = location.href;
    // begin checking for optional info
    if (body.id && !data.action_type) {
        data.action_type = body.id;
    }
    if (body.getAttribute('data-template')) {
        data.download_ad_template = body.getAttribute('data-template');
    }
    if (ads.length) {
        // loop through all ads and grab their zones
        var zones = [];
        for (var i = 0, l = ads.length; i < l; i++) {
            var zone = ads[i].id;
            if (zone) {
                zones.push(zone);
            }
        }
        if (zones.length) {
            data.ad_zones = zones;
        }
    }
    // check the query string for certain data
    if (qs.accel_key) {
        data.ticket = qs.accel_key;
    }
    if (qs.click_id) {
        data.click_id = qs.click_id;
    }
    // include any test data
    if (testString &&
        domainHash &&
        data.active_tests && data.active_tests.length) {
        data.tests = testString;
        data.domain_hash = domainHash;
    }
    // select the correct logging path based on whether this is
    // accelerator or not
    if (qs.accel_key) {
        url += 'accel/';
    }
    // Setup default value for referer
    data.referer = data.referer || doc.referrer;
    // send it all off for data crunching
    $.ajax({
        url: url,
        data: data,
        traditional: true,
        global: false
    });
})();
