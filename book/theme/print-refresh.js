// Remember the chapter that opened the printable book, then return after
// printing or cancelling. Preserve the full URL, including section anchors.
(() => {
    const key = "presolve-print-return";
    const isPrint = (url) => url.pathname.endsWith("/print.html");
    const readReturn = () => {
        try { return window.sessionStorage.getItem(key); }
        catch { return null; }
    };
    const chapterURL = (value) => {
        if (!value) return null;
        try {
            const url = new URL(value, window.location.href);
            return url.origin === window.location.origin && !isPrint(url) ? url.href : null;
        } catch { return null; }
    };

    if (!isPrint(window.location)) {
        document.addEventListener("click", (event) => {
            const link = event.target.closest?.("a[href]");
            if (!link) return;
            const url = new URL(link.href, window.location.href);
            if (url.origin === window.location.origin && isPrint(url)) {
                try { window.sessionStorage.setItem(key, window.location.href); }
                catch { /* Printing remains available without storage. */ }
            }
        }, true);
        return;
    }

    const returnTo = chapterURL(readReturn()) || chapterURL(document.referrer);
    window.addEventListener("afterprint", () => {
        if (!returnTo) return;
        try { window.sessionStorage.removeItem(key); }
        catch { /* The validated referrer can also provide the return URL. */ }
        window.location.replace(returnTo);
    }, { once: true });

    // Direct visits without a known chapter stay on the printable book.
    // Suppress mdBook's one automatic print call on reload in that case too.
    const navigation = window.performance.getEntriesByType("navigation")[0];
    if (!navigation || navigation.type !== "reload") return;
    const print = window.print;
    window.print = function () {
        window.print = print;
    };
})();
