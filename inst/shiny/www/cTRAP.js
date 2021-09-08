/**
 * Clear text selection
 */
function clearSelection() {
    if (window.getSelection) {
        window.getSelection().removeAllRanges();
    } else if (document.selection) {
        document.selection.empty();
    }
}

/**
 * Copy token to clipboard
 */
function copyToken() {
    node = document.getElementById("token");
    if (document.body.createTextRange) {
        const range = document.body.createTextRange();
        range.moveToElementText(node);
        range.select();
        document.execCommand("copy");
        clearSelection();
    } else if (window.getSelection) {
        const selection = window.getSelection();
        const range = document.createRange();
        range.selectNodeContents(node);
        selection.removeAllRanges();
        selection.addRange(range);
        document.execCommand("copy");
        selection.removeAllRanges();
        clearSelection();
    } else {
        console.warn("Could not select text in node: Unsupported browser.");
    }
    return(node)
}
